import numpy as np
import h5py
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge
import json

from typing import Tuple, List

# -----------------------------
# Behavior time alignment + residualization helpers
# -----------------------------

def _require_group_and_dsets(f, group, dsets):
    if group not in f:
        raise KeyError(f"Missing group '{group}' in {f.filename}")
    g = f[group]
    for ds in dsets:
        if ds not in g:
            raise KeyError(f"Missing dataset '{group}/{ds}' in {f.filename}")
    return g


def load_behavior_session_ts(behavior_h5: str, group: str = "dlc_features/session_ts") -> dict:
    """
    Loads session-level behavior time series at 100 Hz (or whatever you stored).
    REQUIRED datasets (as you specified):
      time_s, pos_x_m, pos_y_m, hd_angle_rad, head_ang_vel_rads, body_speed_rms_mps
    Optional:
      hd_valid_mask, pos_valid_mask (if present, used to mask invalid samples)
    """
    import h5py
    with h5py.File(behavior_h5, "r") as f:
        g = _require_group_and_dsets(
            f, group,
            ["time_s", "pos_x_m", "pos_y_m", "hd_angle_rad", "head_ang_vel_rads", "body_speed_rms_mps"]
        )

        out = {
            "time_s": g["time_s"][...].astype(np.float64),
            "pos_x_m": g["pos_x_m"][...].astype(np.float64),
            "pos_y_m": g["pos_y_m"][...].astype(np.float64),
            "hd_angle_rad": g["hd_angle_rad"][...].astype(np.float64),
            "head_ang_vel_rads": g["head_ang_vel_rads"][...].astype(np.float64),
            "body_speed_rms_mps": g["body_speed_rms_mps"][...].astype(np.float64),
        }

        # optional masks
        out["hd_valid_mask"] = g["hd_valid_mask"][...].astype(bool) if "hd_valid_mask" in g else None
        out["pos_valid_mask"] = g["pos_valid_mask"][...].astype(bool) if "pos_valid_mask" in g else None

    # sanity: enforce 1D and same length
    L = len(out["time_s"])
    for k, v in out.items():
        if v is None:
            continue
        if np.asarray(v).ndim != 1:
            raise ValueError(f"{k} must be 1D, got shape {np.asarray(v).shape}")
        if len(v) != L:
            # hard-trim to min length
            L2 = min(L, len(v))
            out["time_s"] = out["time_s"][:L2]
            for kk in ["pos_x_m", "pos_y_m", "hd_angle_rad", "head_ang_vel_rads", "body_speed_rms_mps"]:
                out[kk] = out[kk][:L2]
            if out["hd_valid_mask"] is not None:
                out["hd_valid_mask"] = out["hd_valid_mask"][:L2]
            if out["pos_valid_mask"] is not None:
                out["pos_valid_mask"] = out["pos_valid_mask"][:L2]
            break

    return out


def interp_to_neural_time(neural_t_s: np.ndarray, beh_time_s: np.ndarray, beh_x: np.ndarray) -> np.ndarray:
    """
    Linear interpolation of a 1D behavior vector onto neural times.
    Returns NaNs outside range.
    """
    neural_t_s = np.asarray(neural_t_s, dtype=np.float64)
    beh_time_s = np.asarray(beh_time_s, dtype=np.float64)
    beh_x = np.asarray(beh_x, dtype=np.float64)

    # ensure monotonic behavior time
    order = np.argsort(beh_time_s)
    t = beh_time_s[order]
    x = beh_x[order]

    y = np.interp(neural_t_s, t, x, left=np.nan, right=np.nan)
    return y


def circular_interp_angle(neural_t_s: np.ndarray, beh_time_s: np.ndarray, beh_angle_rad: np.ndarray) -> np.ndarray:
    """
    Interpolate angles by interpolating sin/cos and re-wrapping.
    """
    a = np.asarray(beh_angle_rad, dtype=np.float64)
    s = np.sin(a)
    c = np.cos(a)
    s_i = interp_to_neural_time(neural_t_s, beh_time_s, s)
    c_i = interp_to_neural_time(neural_t_s, beh_time_s, c)
    ang = np.arctan2(s_i, c_i)
    return ang


def build_design_matrix_for_bins(
    t_edges_50ms: np.ndarray | None,
    t_bins: int,
    bin_size_s: float,
    beh: dict,
    *,
    use_terms=("time", "time2", "pos", "hd", "angvel", "speed", "dcenter", "hd_rel_center"),
) -> Tuple[np.ndarray, List[str]]:
    """
    Builds per-bin covariate matrix Xcov for residualization:
      time_s, pos_x, pos_y, hd_angle, head_ang_vel, body_speed, dcenter, hd_rel_center

    use_terms controls which blocks are included.
    """
    # neural bin centers in seconds
    if t_edges_50ms is not None:
        te = np.asarray(t_edges_50ms, dtype=np.float64)
        if te.ndim != 1 or len(te) != t_bins + 1:
            raise ValueError(f"t_edges_50ms must be len(t_bins+1)={t_bins+1}, got {te.shape}")
        t_cent = 0.5 * (te[:-1] + te[1:])
    else:
        t_cent = (np.arange(t_bins, dtype=np.float64) + 0.5) * float(bin_size_s)

    bt = beh["time_s"]

    cols = []
    names = []

    # --- time basis (session drift) ---
    # Normalize time to a stable range so t and t^2 behave well numerically.
    # Maps roughly to [-0.5, 0.5] (clipped to [-1, 1]).
    t0 = np.nanmin(t_cent)
    t1 = np.nanmax(t_cent)
    if np.isfinite(t0) and np.isfinite(t1) and (t1 - t0) > 1e-9:
        t_norm = (t_cent - 0.5 * (t0 + t1)) / (t1 - t0)
        t_norm = np.clip(t_norm, -1.0, 1.0)
    else:
        t_norm = t_cent * 0.0

    if "time" in use_terms:
        cols.append(t_norm.copy())
        names.append("time_norm")

    if "time2" in use_terms:
        cols.append((t_norm ** 2).copy())
        names.append("time_norm2")

    if "pos" in use_terms:
        px = interp_to_neural_time(t_cent, bt, beh["pos_x_m"])
        py = interp_to_neural_time(t_cent, bt, beh["pos_y_m"])
        cols += [px, py]
        names += ["pos_x_m", "pos_y_m"]

    if "hd" in use_terms:
        hd = circular_interp_angle(t_cent, bt, beh["hd_angle_rad"])
        cols.append(hd)
        names.append("hd_angle_rad")

    if "angvel" in use_terms:
        av = interp_to_neural_time(t_cent, bt, beh["head_ang_vel_rads"])
        cols.append(av)
        names.append("head_ang_vel_rads")

    if "speed" in use_terms:
        sp = interp_to_neural_time(t_cent, bt, beh["body_speed_rms_mps"])
        cols.append(sp)
        names.append("body_speed_rms_mps")

    # derived: distance to center
    if "dcenter" in use_terms or "hd_rel_center" in use_terms:
        px = interp_to_neural_time(t_cent, bt, beh["pos_x_m"])
        py = interp_to_neural_time(t_cent, bt, beh["pos_y_m"])
        dcenter = np.sqrt(px**2 + py**2)

    if "dcenter" in use_terms:
        cols.append(dcenter)
        names.append("dist_to_center_m")

    # derived: head direction relative to vector to center
    if "hd_rel_center" in use_terms:
        # angle pointing from animal position to arena center (0,0): atan2(-y, -x)
        ang_to_center = np.arctan2(-py, -px)
        hd = circular_interp_angle(t_cent, bt, beh["hd_angle_rad"])
        # wrap difference to [-pi, pi]
        hd_rel = np.arctan2(np.sin(hd - ang_to_center), np.cos(hd - ang_to_center))
        cols.append(hd_rel)
        names.append("hd_rel_center_rad")

    Xcov = np.column_stack(cols).astype(np.float64)
    return Xcov, names


def residualize_latent_time_series(
    Z_all: np.ndarray,
    Xcov: np.ndarray,
    *,
    ridge_alpha: float = 1.0,
) -> Tuple[np.ndarray, dict]:
    """
    Residualize each latent dimension independently:
      Z_d ~ Xcov  (with intercept)
      Z_res = Z - Z_hat
    Returns (Z_resid, info_dict).
    """
    Z_all = np.asarray(Z_all, dtype=np.float64)
    Xcov = np.asarray(Xcov, dtype=np.float64)
    if Z_all.shape[0] != Xcov.shape[0]:
        L = min(Z_all.shape[0], Xcov.shape[0])
        Z_all = Z_all[:L]
        Xcov = Xcov[:L]

    # valid rows: all covariates finite + Z finite
    m = np.all(np.isfinite(Xcov), axis=1) & np.all(np.isfinite(Z_all), axis=1)
    if m.sum() < 50:
        raise RuntimeError(f"Too few valid rows for residualization: {m.sum()}")

    # z-score covariates on valid rows
    X = Xcov.copy()
    mu = np.nanmean(X[m], axis=0)
    sd = np.nanstd(X[m], axis=0)
    sd[sd == 0] = 1.0
    X = (X - mu) / sd

    model = Ridge(alpha=float(ridge_alpha), fit_intercept=True, random_state=0)
    D = Z_all.shape[1]

    coefs = np.zeros((D, X.shape[1]), dtype=np.float64)
    intercept = np.zeros(D, dtype=np.float64)
    r2 = np.zeros(D, dtype=np.float64)

    Zhat = np.full_like(Z_all, np.nan, dtype=np.float64)

    # fit per-dim
    for d in range(D):
        yd = Z_all[:, d]
        md = m & np.isfinite(yd)
        if md.sum() < 50:
            continue
        model.fit(X[md], yd[md])
        yhat = model.predict(X[md])

        Zhat[md, d] = yhat
        coefs[d] = model.coef_
        intercept[d] = model.intercept_

        # simple R^2 on valid subset
        ss_res = np.sum((yd[md] - yhat) ** 2)
        ss_tot = np.sum((yd[md] - np.mean(yd[md])) ** 2)
        r2[d] = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else np.nan

    Z_res = Z_all - Zhat
    info = {
        "x_mu": mu.astype(np.float32),
        "x_sd": sd.astype(np.float32),
        "coef": coefs.astype(np.float32),
        "intercept": intercept.astype(np.float32),
        "r2_per_dim": r2.astype(np.float32),
        "valid_rows_frac": float(np.mean(m)),
    }
    return Z_res.astype(np.float32), info

# -----------------------------
# Timebase conversion
# -----------------------------
def pulse_periods_to_bin_periods(periods_pulse: np.ndarray,
                                 bins_per_pulse: int = 5,
                                 end_inclusive: bool = True) -> np.ndarray:
    """
    Convert Nx2 pulse periods (start_pulse, end_pulse) to Nx2 50ms-bin periods (start_bin, end_bin).
    Assumes 1 pulse = bins_per_pulse bins.

    If end_inclusive=True: [start_pulse, end_pulse] inclusive.
    If end_inclusive=False: [start_pulse, end_pulse) half-open.

    Returns:
      periods_bin: int array Nx2 with inclusive bin endpoints [start_bin, end_bin]
    """
    periods_pulse = np.asarray(periods_pulse, dtype=np.int64)
    if periods_pulse.ndim != 2 or periods_pulse.shape[1] != 2:
        raise ValueError("periods_pulse must be Nx2")

    start_p = periods_pulse[:, 0]
    end_p = periods_pulse[:, 1]

    start_b = start_p * bins_per_pulse

    if end_inclusive:
        # last pulse contributes bins_per_pulse bins, end bin is end_p*bpp + (bpp-1)
        end_b = end_p * bins_per_pulse + (bins_per_pulse - 1)
    else:
        # end_p is exclusive pulse index, so last included pulse is end_p-1
        end_b = (end_p * bins_per_pulse) - 1

    out = np.stack([start_b, end_b], axis=1)
    return out


def clip_periods(periods_bin: np.ndarray, t_bins: int) -> np.ndarray:
    """Clip periods to [0, t_bins-1] and drop invalid (end < start)."""
    p = np.asarray(periods_bin, dtype=np.int64).copy()
    p[:, 0] = np.clip(p[:, 0], 0, t_bins - 1)
    p[:, 1] = np.clip(p[:, 1], 0, t_bins - 1)
    keep = p[:, 1] >= p[:, 0]
    return p[keep]


def periods_to_mask(periods_bin: np.ndarray, t_bins: int) -> np.ndarray:
    """Create boolean mask length t_bins that is True for bins within any period."""
    mask = np.zeros(t_bins, dtype=bool)
    for s, e in periods_bin:
        mask[s:e+1] = True
    return mask


def extract_fixed_windows_from_periods(periods_bin: np.ndarray,
                                      win_bins: int,
                                      t_bins: int,
                                      align: str = "start") -> np.ndarray:
    """
    For each period [s,e], extract a fixed window of length win_bins aligned to period start or center.

    Returns:
      windows: Nx2 array of bin [w_start, w_end] inclusive of length win_bins
    """
    windows = []
    for s, e in periods_bin:
        if align == "start":
            ws = s
        elif align == "center":
            mid = (s + e) // 2
            ws = mid - win_bins // 2
        else:
            raise ValueError("align must be 'start' or 'center'")
        we = ws + win_bins - 1
        if ws < 0 or we >= t_bins:
            continue
        windows.append((ws, we))
    return np.asarray(windows, dtype=np.int64)


# -----------------------------
# Preprocessing
# -----------------------------
def preprocess_counts(X_counts: np.ndarray,
                      X_is_neurons_by_time: bool = True,
                      transform: str = "sqrt",
                      fit_mask: np.ndarray | None = None):
    """
    Preprocess spike counts:
      - transpose to (t_bins, n_neurons)
      - transform: sqrt or log1p
      - z-score per neuron using bins in fit_mask (recommended: stationary-only)
    Returns:
      Xz: (t_bins, n_neurons) float32
      stats: dict with mu/sd (float32, n_neurons)
    """
    X = X_counts.T if X_is_neurons_by_time else X_counts.copy()
    X = X.astype(np.float64)

    if transform == "sqrt":
        X = np.sqrt(X)
    elif transform == "log1p":
        X = np.log1p(X)
    elif transform is None:
        pass
    else:
        raise ValueError("transform must be 'sqrt', 'log1p', or None")

    t_bins, n_neurons = X.shape

    if fit_mask is None:
        fit_mask = np.ones(t_bins, dtype=bool)
    fit_mask = np.asarray(fit_mask, dtype=bool)
    if fit_mask.sum() < 10:
        raise ValueError("fit_mask has too few bins to compute stable z-scoring")

    mu = X[fit_mask].mean(axis=0, keepdims=True)
    sd = X[fit_mask].std(axis=0, keepdims=True)
    sd[sd == 0] = 1.0
    Xz = (X - mu) / sd

    stats = {"mu": mu.squeeze().astype(np.float32), "sd": sd.squeeze().astype(np.float32)}
    return Xz.astype(np.float32), stats


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

    # residualization (optional)
    behavior_h5: str | None = None,
    do_residualize: bool = False,
    resid_use_terms: tuple = ("time", "pos", "hd", "angvel", "speed", "dcenter", "hd_rel_center"),
    resid_ridge_alpha: float = 1.0,
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

    # Optional: residualize Z_all using session-time behavior covariates
    Z_all_resid = None
    resid_info = None
    if do_residualize:
        if behavior_h5 is None:
            raise ValueError("do_residualize=True requires behavior_h5 path")
        beh = load_behavior_session_ts(behavior_h5, group="dlc_features/session_ts")
        Xcov, xnames = build_design_matrix_for_bins(
            t_edges_50ms=t_edges_50ms,
            t_bins=t_bins,
            bin_size_s=bin_size_s,
            beh=beh,
            use_terms=resid_use_terms,
        )
        Z_all_resid, resid_info = residualize_latent_time_series(
            Z_all, Xcov, ridge_alpha=resid_ridge_alpha
        )
        resid_info["x_names"] = list(xnames)
        resid_info["use_terms"] = list(resid_use_terms)
        resid_info["ridge_alpha"] = float(resid_ridge_alpha)
        resid_info["behavior_h5"] = str(behavior_h5)

    # Extract fixed-length target episode windows
    win_bins = int(round(episode_win_s / bin_size_s))
    tgt_windows = extract_fixed_windows_from_periods(tgt_bin, win_bins=win_bins, t_bins=t_bins, align=align)
    n_ep = tgt_windows.shape[0]
    if n_ep == 0:
        raise RuntimeError("No target windows survived (likely near edges). Adjust episode_win_s or alignment.")

    # Stack all episode trajectories into one array for compact HDF5 storage
    # traj_stack shape: (n_ep * win_bins, n_pc)
    traj_stack = np.empty((n_ep * win_bins, n_pc), dtype=np.float32)
    ep_ptr = np.zeros((n_ep, 2), dtype=np.int64)
    ep_binwin = np.zeros((n_ep, 2), dtype=np.int64)

    traj_stack_resid = None
    if Z_all_resid is not None:
        traj_stack_resid = np.empty((n_ep * win_bins, n_pc), dtype=np.float32)

    row = 0
    k = 0  # number of kept episodes
    kept_orig_idx = []  # maps kept episode index -> original index in tgt_windows

    for i, (bs, be) in enumerate(tgt_windows):
        seg = Z_all[bs:be+1]  # (win_bins, n_pc)
        if seg.shape[0] != win_bins:
            continue
        if not np.all(np.isfinite(seg)):
            continue

        # If residualized trajectories exist, require residual segment to be valid too.
        if traj_stack_resid is not None:
            seg_r = Z_all_resid[bs:be+1]
            if seg_r.shape[0] != win_bins:
                continue
            if not np.all(np.isfinite(seg_r)):
                continue

        # Keep episode: write raw
        traj_stack[row:row+win_bins] = seg

        # Keep episode: write resid (if any)
        if traj_stack_resid is not None:
            traj_stack_resid[row:row+win_bins] = seg_r

        ep_ptr[k] = (row, row + win_bins - 1)
        ep_binwin[k] = (bs, be)
        kept_orig_idx.append(i)

        row += win_bins
        k += 1

    used_ep = int(k)
    traj_stack = traj_stack[:used_ep * win_bins]
    ep_ptr = ep_ptr[:used_ep]
    ep_binwin = ep_binwin[:used_ep]

    if traj_stack_resid is not None:
        traj_stack_resid = traj_stack_resid[:used_ep * win_bins]
        kept_orig_idx = np.asarray(kept_orig_idx, dtype=np.int64)
    else:
        kept_orig_idx = None

    # If any episodes were skipped due to unexpected size, trim
    used_ep = int(row // win_bins)
    traj_stack = traj_stack[:used_ep * win_bins]
    if traj_stack_resid is not None:
        traj_stack_resid = traj_stack_resid[:used_ep * win_bins]

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
        do_residualize=bool(do_residualize),
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

        # Optional residualized latent trace
        if Z_all_resid is not None:
            g_lat.create_dataset("Z_all_resid", data=Z_all_resid, compression="gzip")
            g_res = f.create_group("residualization")
            g_res.attrs["meta_json"] = json.dumps({
                "behavior_h5": resid_info.get("behavior_h5", ""),
                "use_terms": resid_info.get("use_terms", []),
                "x_names": resid_info.get("x_names", []),
                "ridge_alpha": resid_info.get("ridge_alpha", float(resid_ridge_alpha)),
                "valid_rows_frac": resid_info.get("valid_rows_frac", np.nan),
            })
            g_res.create_dataset("x_mu", data=resid_info["x_mu"], compression="gzip")
            g_res.create_dataset("x_sd", data=resid_info["x_sd"], compression="gzip")
            g_res.create_dataset("coef", data=resid_info["coef"], compression="gzip")          # (n_pc, n_cov)
            g_res.create_dataset("intercept", data=resid_info["intercept"], compression="gzip")# (n_pc,)
            g_res.create_dataset("r2_per_dim", data=resid_info["r2_per_dim"], compression="gzip")

        # Episodes
        g_ep = f.create_group("episodes")
        g_ep.create_dataset("target_bin_windows", data=ep_binwin, compression="gzip")     # (n_ep, 2) in 50ms bins
        g_ep.create_dataset("traj_ptr", data=ep_ptr, compression="gzip")                 # (n_ep, 2) row ptr in traj_stack
        g_ep.create_dataset("traj_stack", data=traj_stack, compression="gzip")           # (n_ep*win_bins, n_pc)

        # Optional residualized episode trajectories (same ptr/windowing)
        if traj_stack_resid is not None:
            g_ep.create_dataset("traj_stack_resid", data=traj_stack_resid, compression="gzip")
            # Residualization may drop some episodes; store mapping to original tgt_windows indices.
            g_ep.create_dataset("kept_orig_idx_resid", data=kept_orig_idx, compression="gzip")

    return out_h5_path

