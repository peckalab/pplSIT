import numpy as np
import h5py, json
from typing import Dict, Tuple, Optional


# -----------------------------
# Core metrics
# -----------------------------
def mean_dist_to_centroid(X: np.ndarray) -> float:
    """X: (N,D)"""
    if X.shape[0] < 1:
        return np.nan
    mu = np.mean(X, axis=0, keepdims=True)
    return float(np.mean(np.linalg.norm(X - mu, axis=1)))


def trace_cov(X: np.ndarray) -> float:
    """X: (N,D)"""
    N = X.shape[0]
    if N < 2:
        return np.nan
    Xc = X - np.mean(X, axis=0, keepdims=True)
    cov = (Xc.T @ Xc) / (N - 1)
    return float(np.trace(cov))


def mean_pairwise_distance(X: np.ndarray) -> float:
    """Mean pairwise Euclidean distance for X: (N,D)."""
    N = X.shape[0]
    if N < 2:
        return np.nan
    G = X @ X.T
    sq = np.sum(X**2, axis=1, keepdims=True)
    D2 = sq + sq.T - 2 * G
    D2[D2 < 0] = 0
    iu = np.triu_indices(N, k=1)
    return float(np.mean(np.sqrt(D2[iu])))


def funneling_curve(X_ep: np.ndarray) -> np.ndarray:
    """
    Dispersion(t) curve across episodes:
      disp[t] = mean_e ||x_e,t - mean_e x_e,t||
    X_ep: (N_ep,T,D)
    """
    N, T, D = X_ep.shape
    disp = np.zeros(T, dtype=np.float32)
    for t in range(T):
        Xt = X_ep[:, t, :]
        mu = np.mean(Xt, axis=0, keepdims=True)
        disp[t] = np.mean(np.linalg.norm(Xt - mu, axis=1))
    return disp


def compute_funneling_metrics_from_tensor(
    X_ep: np.ndarray,
    early_bins: np.ndarray,
    late_bins: np.ndarray,
) -> Dict[str, np.ndarray]:
    """
    Computes early/late dispersion metrics + pairwise compression + funneling curve.
    Returns dict with scalars (0-d arrays) + disp_t curve.
    """
    # episode-averaged vectors per window
    X_early = np.mean(X_ep[:, early_bins, :], axis=1)  # (N,D)
    X_late  = np.mean(X_ep[:, late_bins, :], axis=1)   # (N,D)

    # Dispersion A: mean dist to centroid
    dispA_early = mean_dist_to_centroid(X_early)
    dispA_late  = mean_dist_to_centroid(X_late)
    dispA_delta = dispA_late - dispA_early

    # Dispersion B: trace covariance
    dispB_early = trace_cov(X_early)
    dispB_late  = trace_cov(X_late)
    dispB_delta = dispB_late - dispB_early

    # Pairwise distance compression
    pwd_early = mean_pairwise_distance(X_early)
    pwd_late  = mean_pairwise_distance(X_late)
    pwd_delta = pwd_late - pwd_early

    # Funneling curve + index
    disp_t = funneling_curve(X_ep)  # (T,)
    funn_idx = float(np.mean(disp_t[late_bins]) - np.mean(disp_t[early_bins]))

    return {
        "dispA_early": np.array(dispA_early, dtype=np.float64),
        "dispA_late":  np.array(dispA_late,  dtype=np.float64),
        "dispA_delta": np.array(dispA_delta, dtype=np.float64),

        "dispB_early": np.array(dispB_early, dtype=np.float64),
        "dispB_late":  np.array(dispB_late,  dtype=np.float64),
        "dispB_delta": np.array(dispB_delta, dtype=np.float64),

        "pwd_early":   np.array(pwd_early,   dtype=np.float64),
        "pwd_late":    np.array(pwd_late,    dtype=np.float64),
        "pwd_delta":   np.array(pwd_delta,   dtype=np.float64),

        "disp_t":      disp_t.astype(np.float32),
        "funneling_index": np.array(funn_idx, dtype=np.float64),
    }


# -----------------------------
# Nulls
# -----------------------------
def time_shuffle_within_episode(X_ep: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Shuffle time bins within each episode independently."""
    N, T, D = X_ep.shape
    out = np.empty_like(X_ep)
    for e in range(N):
        perm = rng.permutation(T)
        out[e] = X_ep[e, perm, :]
    return out


def build_exclude_mask_from_windows(t_bins: int, windows_bin: np.ndarray) -> np.ndarray:
    """
    windows_bin: (N,2) inclusive [bs,be]
    returns exclude mask length t_bins (True inside any target episode window)
    """
    m = np.zeros(t_bins, dtype=bool)
    for bs, be in windows_bin:
        bs = int(max(0, bs))
        be = int(min(t_bins - 1, be))
        if be >= bs:
            m[bs:be+1] = True
    return m


def sample_random_nontarget_windows(
    Z_all: np.ndarray,
    win_bins: int,
    n_windows: int,
    exclude_mask: np.ndarray,
    rng: np.random.Generator,
    max_tries: int = 200000,
) -> np.ndarray:
    """
    Sample random windows (n_windows, win_bins, D) from Z_all excluding exclude_mask.
    exclude_mask True = forbidden bin.
    """
    T, D = Z_all.shape
    out = np.empty((n_windows, win_bins, D), dtype=np.float32)
    ok = 0
    tries = 0

    while ok < n_windows and tries < max_tries:
        tries += 1
        s = int(rng.integers(0, T - win_bins))
        if np.any(exclude_mask[s:s+win_bins]):
            continue
        out[ok] = Z_all[s:s+win_bins, :]
        ok += 1

    if ok < n_windows:
        raise RuntimeError(f"Could only sample {ok}/{n_windows} random windows. Relax exclusions or increase max_tries.")
    return out


# -----------------------------
# HDF5 IO
# -----------------------------
def _require_del_group(h5: h5py.File, path: str, overwrite: bool) -> h5py.Group:
    if path in h5 and overwrite:
        del h5[path]
    return h5.require_group(path)


def load_episode_tensor_from_core_h5(h5_path: str, use_resid: bool = False) -> Tuple[np.ndarray, np.ndarray, np.ndarray, Dict]:
    """
    Returns:
      X_ep: (n_ep, win_bins, D)
      Z_all: (t_bins, D)
      target_bin_windows: (n_ep,2) inclusive [bs,be]
      meta: dict
    """
    with h5py.File(h5_path, "r") as f:
        meta = json.loads(f.attrs["meta_json"])
        win_bins = int(meta["episode_win_bins"])

        z_key = "latent/Z_all_resid" if use_resid else "latent/Z_all"
        ts_key = "episodes/traj_stack_resid" if use_resid else "episodes/traj_stack"

        if z_key not in f:
            raise KeyError(f"Missing dataset '{z_key}' in {h5_path}. Did you build residualized trajectories?")
        if ts_key not in f:
            raise KeyError(f"Missing dataset '{ts_key}' in {h5_path}. Did you build residualized trajectories?")

        Z_all = f[z_key][...].astype(np.float32)  # (t_bins, D)
        windows = f["episodes/target_bin_windows"][...].astype(np.int64)  # (n_ep,2)
        traj_stack = f[ts_key][...].astype(np.float32)  # (n_ep*win_bins, D)

        n_ep = windows.shape[0]
        D = traj_stack.shape[1]

        # Reconstruct X_ep by reshaping; builder guaranteed fixed-length stacking
        # traj_stack length should be n_ep*win_bins
        if traj_stack.shape[0] != n_ep * win_bins:
            # fall back to using traj_ptr if trimming happened weirdly
            ptr = f["episodes/traj_ptr"][...].astype(np.int64)
            n_ep = ptr.shape[0]
            X_ep = np.empty((n_ep, win_bins, D), dtype=np.float32)
            for i in range(n_ep):
                rs, re = ptr[i]
                seg = traj_stack[rs:re+1]
                if seg.shape[0] != win_bins:
                    raise ValueError("Episode segment length mismatch while reconstructing X_ep from traj_ptr.")
                X_ep[i] = seg
            windows = windows[:n_ep]
        else:
            X_ep = traj_stack.reshape(n_ep, win_bins, D)

    return X_ep, Z_all, windows, meta

def _h5_write_any(g, name, data):
    arr = np.asarray(data)
    if arr.ndim == 0:
        # scalar -> no compression
        g.create_dataset(name, data=arr)
    else:
        g.create_dataset(name, data=arr, compression="gzip")

def save_funneling_results_to_h5(
    out_h5_path: str,
    group: str,
    meta: Dict,
    real: Dict[str, np.ndarray],
    null_time_shuffle: Dict[str, np.ndarray],
    null_random_windows: Dict[str, np.ndarray],
    overwrite_group: bool = True,
):
    with h5py.File(out_h5_path, "a") as f:
        g = _require_del_group(f, group, overwrite_group)
        g.attrs["meta_json"] = json.dumps(meta)

        # real
        greal = g.require_group("real")
        for k, v in real.items():
            _h5_write_any(greal, k, v)

        # shuffle null
        gsh = g.require_group("null_time_shuffle")
        for k, v in null_time_shuffle.items():
            _h5_write_any(gsh, k, v)

        # random windows null
        grd = g.require_group("null_random_windows")
        for k, v in null_random_windows.items():
            _h5_write_any(grd, k, v)


# -----------------------------
# Main per-session runner
# -----------------------------
def compute_and_save_funneling_from_core_h5(
    core_h5_path: str,
    out_h5_path: Optional[str] = None,
    group: str = "funneling",
    use_resid: bool = False,
    # early/late windows in seconds (converted using bin_size_s from meta)
    early_window_s: Tuple[float, float] = (0.0, 1.0),
    late_window_s: Tuple[float, float]  = (4.0, 6.0),
    # null settings
    n_time_shuffles: int = 50,
    n_random_windows_pool: int = 500,   # pool size for random windows
    n_random_bootstrap: int = 50,       # how many bootstrap draws (matched to n_ep)
    seed: int = 0,
    overwrite_group: bool = True,
):
    
    """
    Computes funneling metrics + nulls from a core trajectories H5 produced by build_core_trajectories_h5.
    Saves results to out_h5_path (default: append into core_h5_path).
    """
    if out_h5_path is None:
        out_h5_path = core_h5_path
    
    # Put raw and resid results into separate groups unless user specified otherwise
    if group == "funneling":
        group = "funneling/resid" if use_resid else "funneling/raw"

    rng = np.random.default_rng(int(seed))

    # load core objects
    X_ep, Z_all, tgt_win_bin, meta_in = load_episode_tensor_from_core_h5(core_h5_path, use_resid=use_resid)

    n_ep, win_bins, D = X_ep.shape
    t_bins = Z_all.shape[0]

    bin_size_s = float(meta_in["bin_size_s"])
    # define early/late bins from seconds
    e0 = int(np.floor(early_window_s[0] / bin_size_s))
    e1 = int(np.ceil (early_window_s[1] / bin_size_s))
    l0 = int(np.floor(late_window_s[0]  / bin_size_s))
    l1 = int(np.ceil (late_window_s[1]  / bin_size_s))

    # clip and validate
    e0 = max(0, min(win_bins-1, e0))
    e1 = max(e0+1, min(win_bins, e1))
    l0 = max(0, min(win_bins-1, l0))
    l1 = max(l0+1, min(win_bins, l1))

    early_bins = np.arange(e0, e1, dtype=int)
    late_bins  = np.arange(l0, l1, dtype=int)

    # ---------------- REAL
    real = compute_funneling_metrics_from_tensor(X_ep, early_bins, late_bins)

    # ---------------- NULL 1: within-episode time shuffle
    shuf_vals = {k: [] for k in real.keys() if k != "disp_t"}
    shuf_disp_t = []

    for _ in range(int(n_time_shuffles)):
        Xs = time_shuffle_within_episode(X_ep, rng)
        m = compute_funneling_metrics_from_tensor(Xs, early_bins, late_bins)
        for k, v in m.items():
            if k == "disp_t":
                shuf_disp_t.append(v)
            else:
                shuf_vals[k].append(float(v))

    null_time_shuffle = {
        "disp_t_values": np.asarray(shuf_disp_t, dtype=np.float32),  # (n_shuf, T)
    }
    for k, vals in shuf_vals.items():
        vals = np.asarray(vals, dtype=np.float64)
        null_time_shuffle[f"{k}_values"] = vals
        null_time_shuffle[f"{k}_mean"] = np.array(np.mean(vals), dtype=np.float64)
        null_time_shuffle[f"{k}_std"]  = np.array(np.std(vals, ddof=1) if len(vals) > 1 else np.nan, dtype=np.float64)

    # ---------------- NULL 2: random non-target windows
    # Exclude any bins that are part of target windows
    exclude_mask = build_exclude_mask_from_windows(t_bins, tgt_win_bin)

    # Pool of random windows
    pool = sample_random_nontarget_windows(
        Z_all=Z_all,
        win_bins=win_bins,
        n_windows=int(n_random_windows_pool),
        exclude_mask=exclude_mask,
        rng=rng,
    )  # (pool, win_bins, D)

    rand_vals = {k: [] for k in real.keys() if k != "disp_t"}
    rand_disp_t = []

    for _ in range(int(n_random_bootstrap)):
        idx = rng.choice(pool.shape[0], size=n_ep, replace=True)
        Xr = pool[idx]  # (n_ep, win_bins, D)
        m = compute_funneling_metrics_from_tensor(Xr, early_bins, late_bins)
        for k, v in m.items():
            if k == "disp_t":
                rand_disp_t.append(v)
            else:
                rand_vals[k].append(float(v))

    null_random_windows = {
        "disp_t_values": np.asarray(rand_disp_t, dtype=np.float32),  # (n_boot, T)
    }
    for k, vals in rand_vals.items():
        vals = np.asarray(vals, dtype=np.float64)
        null_random_windows[f"{k}_values"] = vals
        null_random_windows[f"{k}_mean"] = np.array(np.mean(vals), dtype=np.float64)
        null_random_windows[f"{k}_std"]  = np.array(np.std(vals, ddof=1) if len(vals) > 1 else np.nan, dtype=np.float64)

    # ---------------- SAVE
    meta = {
        "source_core_h5": str(core_h5_path),
        "bin_size_s": float(bin_size_s),
        "n_ep": int(n_ep),
        "win_bins": int(win_bins),
        "D": int(D),
        "early_window_s": list(map(float, early_window_s)),
        "late_window_s": list(map(float, late_window_s)),
        "early_bins": [int(e0), int(e1)],
        "late_bins": [int(l0), int(l1)],
        "n_time_shuffles": int(n_time_shuffles),
        "n_random_windows_pool": int(n_random_windows_pool),
        "n_random_bootstrap": int(n_random_bootstrap),
        "seed": int(seed),
        "use_resid": bool(use_resid),
        "note": "Real computed on target episodes; random windows sampled from latent/Z_all excluding episodes/target_bin_windows.",
    }

    save_funneling_results_to_h5(
        out_h5_path=out_h5_path,
        group=group,
        meta=meta,
        real=real,
        null_time_shuffle=null_time_shuffle,
        null_random_windows=null_random_windows,
        overwrite_group=overwrite_group,
    )

    return out_h5_path
