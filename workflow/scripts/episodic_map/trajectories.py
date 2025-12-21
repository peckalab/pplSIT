import os, sys
import h5py
import json
import numpy as np
from scipy.ndimage import gaussian_filter1d


# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.episodic_map.trajectories import *

# -----------------------------
# Main builder
# -----------------------------
def build_core_trajectories_h5(
    out_h5_path: str,
    X_counts_50ms: np.ndarray,
    state_periods_pulse: dict,
    t_edges_50ms: np.ndarray | None = None,
    *,
    X_is_neurons_by_time: bool = True,
    bins_per_pulse: int = 5,           # 250ms / 50ms = 5
    pulse_end_inclusive: bool = True,

    # which states to use
    key_target: str  = "tgt_sta_succ_mx",
    key_bgr_sta: str = "bgr_sta_mx",
    key_sil_sta: str = "sil_sta_mx",

    # PCA fit choices
    transform: str = "sqrt",
    pca_n_components: int = 20,

    # episode windowing
    episode_win_s: float = 6.0,       # 6 s target episodes
    bin_size_s: float = 0.05,         # 50 ms
    align: str = "start",             # align to target period start
):
    """
    Builds and stores:
      - stationary-fit z-scoring stats
      - PCA fit on stationary-only bins (bgr_sta + sil_sta, excluding target)
      - projections for all bins
      - extracted fixed-length target episode trajectories in PCA space

    state_periods_pulse: dict of Nx2 arrays in PULSE INDICES (4Hz), each row (start_pulse, end_pulse).
    """
    # Prepare X in (t_bins, n_neurons)
    X = X_counts_50ms.T if X_is_neurons_by_time else X_counts_50ms
    X = np.asarray(X)
    t_bins, n_neurons = X.shape

    # Convert state periods from pulse->bin
    def get_periods_bin(key):
        if key not in state_periods_pulse:
            raise KeyError(f"Missing state periods key: {key}")
        p_pulse = np.asarray(state_periods_pulse[key], dtype=np.int64)
        p_bin = pulse_periods_to_bin_periods(
            p_pulse, bins_per_pulse=bins_per_pulse, end_inclusive=pulse_end_inclusive
        )
        p_bin = clip_periods(p_bin, t_bins=t_bins)
        return p_bin

    tgt_bin = get_periods_bin(key_target)
    bgr_bin = get_periods_bin(key_bgr_sta)
    sil_bin = get_periods_bin(key_sil_sta)

    # Build masks
    mask_tgt = periods_to_mask(tgt_bin, t_bins)
    mask_sta = periods_to_mask(bgr_bin, t_bins) | periods_to_mask(sil_bin, t_bins)

    # Exclude target from stationary fit (important)
    mask_pca_fit = mask_sta & (~mask_tgt)

    # Preprocess + z-score using stationary-fit bins
    Xz, zstats = preprocess_counts(
        X_counts_50ms, X_is_neurons_by_time=X_is_neurons_by_time,
        transform=transform, fit_mask=mask_pca_fit
    )

    # PCA fit on stationary-only (excluding target)
    n_pc = int(min(pca_n_components, n_neurons))
    pca = PCA(n_components=n_pc, random_state=0)
    pca.fit(Xz[mask_pca_fit])

    # Project all bins
    Z_all = pca.transform(Xz).astype(np.float32)  # (t_bins, n_pc)

    # Extract fixed-length target episode windows
    win_bins = int(round(episode_win_s / bin_size_s))
    tgt_windows = extract_fixed_windows_from_periods(tgt_bin, win_bins=win_bins, t_bins=t_bins, align=align)
    n_ep = tgt_windows.shape[0]
    if n_ep == 0:
        raise RuntimeError("No target windows survived (likely near edges). Adjust episode_win_s or alignment.")

    # Stack all episode trajectories into one array for compact HDF5 storage
    # traj_stack shape: (n_ep * win_bins, n_pc)
    traj_stack = np.empty((n_ep * win_bins, n_pc), dtype=np.float32)
    ep_ptr = np.zeros((n_ep, 2), dtype=np.int64)  # start_row, end_row inclusive in traj_stack
    ep_binwin = np.zeros((n_ep, 2), dtype=np.int64)  # start_bin, end_bin inclusive

    row = 0
    for i, (bs, be) in enumerate(tgt_windows):
        seg = Z_all[bs:be+1]  # (win_bins, n_pc)
        if seg.shape[0] != win_bins:
            continue
        traj_stack[row:row+win_bins] = seg
        ep_ptr[i] = (row, row + win_bins - 1)
        ep_binwin[i] = (bs, be)
        row += win_bins

    # If any episodes were skipped due to unexpected size, trim
    used_ep = int(row // win_bins)
    traj_stack = traj_stack[:used_ep * win_bins]
    ep_ptr = ep_ptr[:used_ep]
    ep_binwin = ep_binwin[:used_ep]

    # Write HDF5
    meta = dict(
        bin_size_s=float(bin_size_s),
        bins_per_pulse=int(bins_per_pulse),
        pulse_end_inclusive=bool(pulse_end_inclusive),
        transform=transform,
        pca_n_components=int(n_pc),
        episode_win_s=float(episode_win_s),
        episode_win_bins=int(win_bins),
        align=align,
        keys=dict(target=key_target, bgr_sta=key_bgr_sta, sil_sta=key_sil_sta),
        n_neurons=int(n_neurons),
        t_bins=int(t_bins),
        n_target_episodes=int(used_ep),
    )

    with h5py.File(out_h5_path, "w") as f:
        f.attrs["meta_json"] = json.dumps(meta)

        g_in = f.create_group("inputs")
        if t_edges_50ms is not None:
            g_in.create_dataset("t_edges_50ms", data=np.asarray(t_edges_50ms, dtype=np.float64), compression="gzip")

        # Store state periods in bin units
        g_state = f.create_group("states")
        g_state.create_dataset("target_periods_bin", data=tgt_bin, compression="gzip")
        g_state.create_dataset("bgr_sta_periods_bin", data=bgr_bin, compression="gzip")
        g_state.create_dataset("sil_sta_periods_bin", data=sil_bin, compression="gzip")
        g_state.create_dataset("mask_tgt", data=mask_tgt.astype(np.uint8), compression="gzip")
        g_state.create_dataset("mask_sta", data=mask_sta.astype(np.uint8), compression="gzip")
        g_state.create_dataset("mask_pca_fit", data=mask_pca_fit.astype(np.uint8), compression="gzip")

        # Preprocessing stats
        g_pp = f.create_group("preproc")
        g_pp.create_dataset("z_mu", data=zstats["mu"], compression="gzip")
        g_pp.create_dataset("z_sd", data=zstats["sd"], compression="gzip")

        # PCA model
        g_pca = f.create_group("pca")
        g_pca.create_dataset("components", data=pca.components_.astype(np.float32), compression="gzip")  # (n_pc, n_neurons)
        g_pca.create_dataset("mean", data=pca.mean_.astype(np.float32), compression="gzip")
        g_pca.create_dataset("explained_variance", data=pca.explained_variance_.astype(np.float32), compression="gzip")
        g_pca.create_dataset("explained_variance_ratio", data=pca.explained_variance_ratio_.astype(np.float32), compression="gzip")

        # Projections (all time bins)
        g_lat = f.create_group("latent")
        g_lat.create_dataset("Z_all", data=Z_all, compression="gzip")  # (t_bins, n_pc)

        # Episodes
        g_ep = f.create_group("episodes")
        g_ep.create_dataset("target_bin_windows", data=ep_binwin, compression="gzip")     # (n_ep, 2) in 50ms bins
        g_ep.create_dataset("traj_ptr", data=ep_ptr, compression="gzip")                 # (n_ep, 2) row ptr in traj_stack
        g_ep.create_dataset("traj_stack", data=traj_stack, compression="gzip")           # (n_ep*win_bins, n_pc)

    return out_h5_path


cfg = snakemake.config['trajectories']

STATE_KEYS = [
    "tgt_sta_succ_mx",
    "bgr_sta_mx",
    "sil_sta_mx",
]

# -----------------------------
# Reading datasets
# -----------------------------
meta_file = snakemake.input[0]
segm_file = snakemake.input[1]
actm_file = snakemake.input[2]

with h5py.File(meta_file, 'r') as f:
    tl = np.array(f['processed']['timeline'])
    tgt_mx = np.array(f['processed']['target_matrix'])
    sound_events = np.array(f['processed']['sound_events'])

with h5py.File(actm_file, 'r') as f:
    X_counts_50ms = np.array(f['mx_50ms']['mx'])
    t_edges_50ms  = np.array(f['mx_50ms']['bins'])
    
#X_counts_50ms = gaussian_filter1d(X_counts_50ms, sigma=5, axis=0, mode="nearest")

with h5py.File(segm_file, 'r') as f:
    state_periods_pulse = {}
    for idxs_name in STATE_KEYS:
        state_periods_pulse[idxs_name] = np.array(f[idxs_name])[:, :2].astype(np.int32)

out = build_core_trajectories_h5(
    out_h5_path=snakemake.output[0],
    X_counts_50ms=X_counts_50ms,
    state_periods_pulse=state_periods_pulse,
    t_edges_50ms=t_edges_50ms,
    X_is_neurons_by_time=True,      # set correctly for your matrix
    episode_win_s=cfg['episode_win_s'],
    pca_n_components=cfg['pca_n_components'],
    transform=cfg['transform'],
)




# -----------------------------
# Helper to read back a single episode trajectory
# -----------------------------
# def load_episode_traj(h5_path: str, ep_idx: int):
#     with h5py.File(h5_path, "r") as f:
#         ptr = f["episodes/traj_ptr"][ep_idx]
#         Z = f["episodes/traj_stack"][ptr[0]:ptr[1]+1]
#     return Z
