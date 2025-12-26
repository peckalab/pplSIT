import numpy as np
import h5py, json
from typing import Dict, Tuple, Optional, Callable


# ============================================================
# Helpers: IO (trajectories.h5)
# ============================================================
def load_core_episodes_latent(trajectories_h5: str):
    """
    Loads latent episodes from your build_core_trajectories_h5 output.
    Returns:
      Z_all: (t_bins, D)
      target_bin_windows: (n_ep, 2) inclusive [bs, be] in 50ms bins
      X_ep_lat: (n_ep, T, D) episode tensor in latent space
      meta: dict (from attrs["meta_json"])
    """
    with h5py.File(trajectories_h5, "r") as f:
        meta = json.loads(f.attrs["meta_json"])
        win_bins = int(meta["episode_win_bins"])

        Z_all = f["latent/Z_all"][...].astype(np.float32)
        target_bin_windows = f["episodes/target_bin_windows"][...].astype(np.int64)  # inclusive
        traj_stack = f["episodes/traj_stack"][...].astype(np.float32)

        n_ep = target_bin_windows.shape[0]
        D = traj_stack.shape[1]

        # Reconstruct tensor (the builder stacks fixed windows in order)
        if traj_stack.shape[0] == n_ep * win_bins:
            X_ep_lat = traj_stack.reshape(n_ep, win_bins, D)
        else:
            # fallback using traj_ptr
            ptr = f["episodes/traj_ptr"][...].astype(np.int64)  # inclusive row ptrs
            n_ep2 = ptr.shape[0]
            X_ep_lat = np.empty((n_ep2, win_bins, D), dtype=np.float32)
            for i in range(n_ep2):
                rs, re = ptr[i]
                seg = traj_stack[rs:re+1]
                if seg.shape[0] != win_bins:
                    raise ValueError("traj_ptr segment length mismatch; episode windows not fixed-length.")
                X_ep_lat[i] = seg
            target_bin_windows = target_bin_windows[:n_ep2]
            n_ep = n_ep2

    return Z_all, target_bin_windows, X_ep_lat, meta


# ============================================================
# Helpers: IO (behavior_features.h5)
# ============================================================
import numpy as np
import h5py

def load_episode_behavior_features(
    beh_h5: str,
    *,
    require_position: bool = True,
    targets_periods_times: np.ndarray | None = None,
    pos_valid_thr: float = 0.0,  # if you ever store fractional validity, else ignored
):
    """
    Loads behavior features from a behavior_features.h5 with:
      - dlc_features/session_ts group containing:
        ['body_speed_rms_mps', 'hd_angle_rad', 'hd_valid_mask', 'head_ang_vel_rads',
         'pos_valid_mask', 'pos_x_m', 'pos_y_m', 'time_s']
      - dlc_features/episode_feats group containing:
        ['hd_mean_rad', 'hd_valid_frac', 'hd_var', 'head_turn_rate_rads', 'stillness_mps']

    Additionally:
      - Requires position to be present (if require_position=True).
      - Computes episode mean position (pos_x, pos_y) by averaging session_ts pos within each episode
        using targets_periods_times (N_ep x 2, seconds). If targets_periods_times is None, tries to
        load it from dlc_features/episode_feats/targets_periods_times if present.

    Returns dict epf with keys:
      hd_mean, hd_var, still, turn, hd_valid,
      pos_x, pos_y,
      and optionally ep_times.
    """

    epf = {}

    with h5py.File(beh_h5, "r") as f:
        # --- required groups
        if "dlc_features/session_ts" not in f:
            raise KeyError(f"Missing group 'dlc_features/session_ts' in {beh_h5}")
        if "dlc_features/episode_feats" not in f:
            raise KeyError(f"Missing group 'dlc_features/episode_feats' in {beh_h5}")

        g_ts = f["dlc_features/session_ts"]
        g_ep = f["dlc_features/episode_feats"]

        # --- episode-level direct features (required)
        req_ep = ["hd_mean_rad", "hd_valid_frac", "hd_var", "head_turn_rate_rads", "stillness_mps"]
        missing_ep = [k for k in req_ep if k not in g_ep]
        if missing_ep:
            raise KeyError(f"Missing episode_feats datasets: {missing_ep}")

        epf["hd_mean"]  = np.asarray(g_ep["hd_mean_rad"][...])
        epf["hd_valid"] = np.asarray(g_ep["hd_valid_frac"][...])
        epf["hd_var"]   = np.asarray(g_ep["hd_var"][...])
        epf["turn"]     = np.asarray(g_ep["head_turn_rate_rads"][...])
        epf["still"]    = np.asarray(g_ep["stillness_mps"][...])
        
        # NEW: center distance and center-referenced HD (episode-level)
        if "dist_to_center_m" in g_ep:
            epf["dist_center"] = np.asarray(g_ep["dist_to_center_m"][...])

        # depending on your naming in behavior_features.h5
        if "hd_rel_center_mean_rad" in g_ep:
            epf["hd_rel_center_mean"] = np.asarray(g_ep["hd_rel_center_mean_rad"][...])
        elif "hd_rel_center_mean_rad" in g_ep:
            epf["hd_rel_center_mean"] = np.asarray(g_ep["hd_rel_center_mean_rad"][...])

        # --- episode times: provided or load if present
        if targets_periods_times is None:
            if "targets_periods_times" in g_ep:
                targets_periods_times = np.asarray(g_ep["targets_periods_times"][...])
        if targets_periods_times is not None:
            targets_periods_times = np.asarray(targets_periods_times, dtype=float)
            if targets_periods_times.ndim != 2 or targets_periods_times.shape[1] != 2:
                raise ValueError("targets_periods_times must be (n_ep, 2) [t_start, t_end] in seconds.")
            epf["ep_times"] = targets_periods_times

        # --- session time series position (required if require_position)
        req_ts = ["pos_x_m", "pos_y_m", "time_s", "pos_valid_mask"]
        missing_ts = [k for k in req_ts if k not in g_ts]
        if missing_ts and require_position:
            raise KeyError(f"Missing session_ts datasets required for position: {missing_ts}")

        if require_position:
            t_s = np.asarray(g_ts["time_s"][...], dtype=float)
            pos_x = np.asarray(g_ts["pos_x_m"][...], dtype=float)
            pos_y = np.asarray(g_ts["pos_y_m"][...], dtype=float)
            pos_valid = np.asarray(g_ts["pos_valid_mask"][...]).astype(bool)

            # basic alignment guard: clip all to same length
            L = min(len(t_s), len(pos_x), len(pos_y), len(pos_valid))
            t_s = t_s[:L]; pos_x = pos_x[:L]; pos_y = pos_y[:L]; pos_valid = pos_valid[:L]

            # require episode times to compute episode mean position
            if targets_periods_times is None:
                raise ValueError(
                    "Position is required, but targets_periods_times was not provided and not found in file. "
                    "Pass targets_periods_times (n_ep,2) seconds."
                )

            n_ep = targets_periods_times.shape[0]
            if epf["hd_mean"].shape[0] != n_ep:
                # allow trimming to min as a safe fallback
                n = min(n_ep, epf["hd_mean"].shape[0])
                targets_periods_times = targets_periods_times[:n]
                for k in ["hd_mean","hd_valid","hd_var","turn","still"]:
                    epf[k] = epf[k][:n]
                n_ep = n

            pos_x_ep = np.full(n_ep, np.nan, dtype=float)
            pos_y_ep = np.full(n_ep, np.nan, dtype=float)
            pos_valid_frac_ep = np.zeros(n_ep, dtype=float)

            # fast indexing via searchsorted on monotonic time
            # assumes time_s is sorted
            for i in range(n_ep):
                t0, t1 = targets_periods_times[i]
                if not np.isfinite(t0) or not np.isfinite(t1) or (t1 <= t0):
                    continue
                a = int(np.searchsorted(t_s, t0, side="left"))
                b = int(np.searchsorted(t_s, t1, side="right"))
                a = max(0, min(a, L))
                b = max(0, min(b, L))
                if b - a < 2:
                    continue

                m = pos_valid[a:b]
                pos_valid_frac_ep[i] = float(np.mean(m)) if (b > a) else 0.0
                if pos_valid_frac_ep[i] <= pos_valid_thr:
                    continue

                xx = pos_x[a:b].copy()
                yy = pos_y[a:b].copy()
                xx[~m] = np.nan
                yy[~m] = np.nan

                pos_x_ep[i] = float(np.nanmean(xx))
                pos_y_ep[i] = float(np.nanmean(yy))

            # Require that we actually computed position for most episodes
            if np.sum(np.isfinite(pos_x_ep) & np.isfinite(pos_y_ep)) < max(3, int(0.3 * n_ep)):
                raise RuntimeError(
                    "Episode mean position could not be computed for enough episodes. "
                    "Check time alignment, pos_valid_mask, or targets_periods_times."
                )

            epf["pos_x"] = pos_x_ep.astype(np.float32)
            epf["pos_y"] = pos_y_ep.astype(np.float32)
            epf["pos_valid_frac"] = pos_valid_frac_ep.astype(np.float32)

    return epf


# ============================================================
# Math: similarity metrics (modular)
# ============================================================
def _normalize_rows(X: np.ndarray, eps: float = 1e-9) -> np.ndarray:
    """X: (..., D) normalize last axis to unit norm."""
    n = np.linalg.norm(X, axis=-1, keepdims=True)
    return X / (n + eps)


def cosine_similarity_phase_aligned(X_ep: np.ndarray, time_bins: np.ndarray) -> np.ndarray:
    """
    Phase-aligned cosine similarity averaged over selected time bins.
    X_ep: (n_ep, T, D)
    time_bins: 1D indices into T
    Returns S: (n_ep, n_ep), symmetric, diag=1.
    """
    X = X_ep[:, time_bins, :]  # (n_ep, Tw, D)
    X = _normalize_rows(X)     # normalize each vector

    n_ep, Tw, D = X.shape
    S = np.zeros((n_ep, n_ep), dtype=np.float32)
    for t in range(Tw):
        Xt = X[:, t, :]  # (n_ep, D)
        S += (Xt @ Xt.T).astype(np.float32)
    S /= float(Tw)
    # numerical tidy
    S = np.clip(S, -1.0, 1.0)
    np.fill_diagonal(S, 1.0)
    return S


def neg_euclidean_similarity_phase_aligned(X_ep: np.ndarray, time_bins: np.ndarray) -> np.ndarray:
    """
    Phase-aligned similarity defined as negative Euclidean distance (averaged over time).
    Similarity more positive => closer.
    """
    X = X_ep[:, time_bins, :]  # (n_ep, Tw, D)
    n_ep, Tw, D = X.shape
    S = np.zeros((n_ep, n_ep), dtype=np.float32)
    for t in range(Tw):
        Xt = X[:, t, :]  # (n_ep, D)
        # pairwise distances via (x-y)^2 = x^2 + y^2 - 2x·y
        G = Xt @ Xt.T
        sq = np.sum(Xt**2, axis=1, keepdims=True)
        D2 = sq + sq.T - 2*G
        D2[D2 < 0] = 0
        S += (-np.sqrt(D2)).astype(np.float32)
    S /= float(Tw)
    np.fill_diagonal(S, 0.0)  # distance is 0 => similarity 0 on diag
    return S


# ============================================================
# PV episode tensor extraction (Option B)
# ============================================================
def extract_episode_tensor_from_counts(
    X_counts_50ms: np.ndarray,
    target_bin_windows: np.ndarray,
    win_bins: int,
    X_is_neurons_by_time: bool = True,
    transform: str = "sqrt",
    zscore: bool = True,
    zfit_mask: Optional[np.ndarray] = None,
    eps: float = 1e-9,
) -> np.ndarray:
    """
    Build (n_ep, win_bins, n_neurons) tensor from raw counts at 50ms bins.
    Uses the windows from trajectories.h5 (inclusive).
    Preprocessing options are modular (sqrt + zscore).

    Notes:
      - zfit_mask if provided should be length t_bins boolean mask for fitting mu/sd.
        If None: fit on all bins.
    """
    X = X_counts_50ms.T if X_is_neurons_by_time else X_counts_50ms
    X = np.asarray(X, float)               # (t_bins, n_neurons)
    t_bins, n_neu = X.shape

    # transform
    if transform == "sqrt":
        Xp = np.sqrt(np.maximum(X, 0.0))
    elif transform == "log1p":
        Xp = np.log1p(np.maximum(X, 0.0))
    elif transform in (None, "none"):
        Xp = X
    else:
        raise ValueError(f"Unknown transform: {transform}")

    # z-score across time (fit on mask)
    if zscore:
        if zfit_mask is None:
            zfit_mask = np.ones(t_bins, dtype=bool)
        mu = np.mean(Xp[zfit_mask], axis=0)
        sd = np.std(Xp[zfit_mask], axis=0, ddof=1)
        sd = np.where(sd < eps, 1.0, sd)
        Xp = (Xp - mu) / sd

    # episodes tensor
    n_ep = target_bin_windows.shape[0]
    X_ep = np.empty((n_ep, win_bins, n_neu), dtype=np.float32)
    keep = np.ones(n_ep, dtype=bool)

    for i in range(n_ep):
        bs, be = target_bin_windows[i]
        if be - bs + 1 != win_bins:
            keep[i] = False
            continue
        if bs < 0 or be >= t_bins:
            keep[i] = False
            continue
        X_ep[i] = Xp[bs:be+1].astype(np.float32)

    if not np.all(keep):
        X_ep = X_ep[keep]
    return X_ep


# ============================================================
# Drivers: pairwise covariate matrices
# ============================================================
def circ_dist(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Circular distance between angles in radians; returns in [0, pi]."""
    d = np.angle(np.exp(1j*(a[:, None] - b[None, :])))
    return np.abs(d).astype(np.float32)


def pairwise_absdiff(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, float)
    return np.abs(x[:, None] - x[None, :]).astype(np.float32)


def pairwise_euclid_xy(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    x = np.asarray(x, float); y = np.asarray(y, float)
    dx = x[:, None] - x[None, :]
    dy = y[:, None] - y[None, :]
    return np.sqrt(dx*dx + dy*dy).astype(np.float32)


def pairwise_time_sep_from_windows(target_bin_windows: np.ndarray, bin_size_s: float) -> np.ndarray:
    """
    Use episode midpoints in bin units to get |Δt| in seconds.
    windows inclusive [bs,be].
    """
    mid = 0.5*(target_bin_windows[:, 0] + target_bin_windows[:, 1]).astype(float)
    mid_s = mid * float(bin_size_s)
    return pairwise_absdiff(mid_s)


# ============================================================
# Vectorization for partial corr later
# ============================================================
def upper_tri_vector(M: np.ndarray) -> np.ndarray:
    """Vectorize upper triangle (k=1)."""
    iu = np.triu_indices(M.shape[0], k=1)
    return M[iu]


def upper_tri_mask_from_valid(valid: np.ndarray) -> np.ndarray:
    """
    valid: (n_ep,) boolean indicating whether episode has valid features.
    Returns mask over upper-tri vectorization.
    """
    n = len(valid)
    iu = np.triu_indices(n, k=1)
    return (valid[iu[0]] & valid[iu[1]])


# ============================================================
# Main pipeline: compute + save
# ============================================================
def compute_episode_similarity_and_drivers_to_h5(
    out_h5_path: str,
    trajectories_h5: str,
    behavior_h5: str,
    # optional PV counts (for option B)
    X_counts_50ms: Optional[np.ndarray] = None,
    X_is_neurons_by_time: bool = True,
    # episode selection (indices into trajectories episodes)
    ep_select_idx: Optional[np.ndarray] = None,
    # similarity choices
    latent_similarity: str = "cosine",  # "cosine" or "neg_euclid"
    pv_similarity: str = "cosine",      # "cosine" (recommended)
    # phase window
    strip_first_bins: int = 10,         # e.g. 0.5 s at 50 ms
    use_time_bins: Optional[np.ndarray] = None,  # overrides computed time bins
    # PV preprocessing
    pv_transform: str = "sqrt",
    pv_zscore: bool = True,
    # saving
    group: str = "episode_similarity",
    overwrite: bool = True,
):
    """
    Produces and saves:
      - similarity matrices: latent + PV
      - covariate matrices: space/time/HD/stillness/turn
      - vectorized upper-tri arrays + masks for partial-corr later
      - meta json

    Requires:
      trajectories_h5: core trajectories file
      behavior_h5: behavior features file
      X_counts_50ms: raw spike counts at 50ms (optional; if None, PV similarity skipped)
    """

    # ---- Load latent episodes + windows
    Z_all, target_bin_windows, X_ep_lat, meta_traj = load_core_episodes_latent(trajectories_h5)
    win_bins = int(meta_traj["episode_win_bins"])
    bin_size_s = float(meta_traj["bin_size_s"])

    n_ep_all = target_bin_windows.shape[0]

    # ---- Episode selection
    if ep_select_idx is None:
        ep_select_idx = np.arange(n_ep_all, dtype=int)
    else:
        ep_select_idx = np.asarray(ep_select_idx, dtype=int)

    target_bin_windows = target_bin_windows[ep_select_idx]
    X_ep_lat = X_ep_lat[ep_select_idx]
    n_ep = X_ep_lat.shape[0]

    # ---- Phase bins
    if use_time_bins is None:
        tb = np.arange(win_bins, dtype=int)
        tb = tb[tb >= int(strip_first_bins)]
        if tb.size == 0:
            raise ValueError("No time bins left after strip_first_bins.")
        use_time_bins = tb
    else:
        use_time_bins = np.asarray(use_time_bins, dtype=int)

    # ---- Latent similarity
    if latent_similarity == "cosine":
        S_lat = cosine_similarity_phase_aligned(X_ep_lat, use_time_bins)
    elif latent_similarity == "neg_euclid":
        S_lat = neg_euclidean_similarity_phase_aligned(X_ep_lat, use_time_bins)
    else:
        raise ValueError(f"Unknown latent_similarity: {latent_similarity}")

    # ---- PV similarity (optional)
    S_pv = None
    if X_counts_50ms is not None:
        # Fit z-score on stationary mask if available in trajectories.h5, else all bins
        zfit_mask = None
        with h5py.File(trajectories_h5, "r") as f:
            if "states/mask_pca_fit" in f:
                zfit_mask = f["states/mask_pca_fit"][...].astype(bool)

        X_ep_pv = extract_episode_tensor_from_counts(
            X_counts_50ms=X_counts_50ms,
            target_bin_windows=target_bin_windows,
            win_bins=win_bins,
            X_is_neurons_by_time=X_is_neurons_by_time,
            transform=pv_transform,
            zscore=pv_zscore,
            zfit_mask=zfit_mask,
        )

        # Important: if any episodes were dropped due to window mismatch, align by truncation
        n_ep2 = X_ep_pv.shape[0]
        if n_ep2 != n_ep:
            # safest: keep only first n_ep2 episodes (or you can pass a consistent ep_select_idx)
            n_keep = min(n_ep, n_ep2)
            X_ep_lat = X_ep_lat[:n_keep]
            target_bin_windows = target_bin_windows[:n_keep]
            S_lat = S_lat[:n_keep, :n_keep]
            X_ep_pv = X_ep_pv[:n_keep]
            n_ep = n_keep

        if pv_similarity == "cosine":
            S_pv = cosine_similarity_phase_aligned(X_ep_pv, use_time_bins)
        else:
            raise ValueError(f"Unknown pv_similarity: {pv_similarity}")

    # ---- Load episode-level behavior features
    epf = load_episode_behavior_features(behavior_h5)

    # Required features (HD/stillness/turn)
    hd_mean = epf.get("hd_mean", None)
    still = epf.get("still", None)
    turn = epf.get("turn", None)
    hd_valid = epf.get("hd_valid", None)

    # Optional episode mean position
    pos_x = epf.get("pos_x", None)
    pos_y = epf.get("pos_y", None)
    dist_center = epf.get("dist_center", None)
    hd_rel_center_mean = epf.get("hd_rel_center_mean", None)

    # Align behavior episode features with selected episodes
    # (Assumes behavior ep ordering matches trajectories ep ordering; if not, you will pass ep_select_idx accordingly.)
    def sel(v):
        if v is None:
            return None
        v = np.asarray(v)
        return v[ep_select_idx][:n_ep]  # also handle truncation

    hd_mean = sel(hd_mean)
    still = sel(still)
    turn = sel(turn)
    hd_valid = sel(hd_valid)
    pos_x = sel(pos_x)
    pos_y = sel(pos_y)

    # ---- Covariate matrices
    cov = {}

    cov["dt_s"] = pairwise_time_sep_from_windows(target_bin_windows, bin_size_s=bin_size_s)

    if hd_mean is not None:
        cov["dhd_rad"] = circ_dist(hd_mean.astype(float), hd_mean.astype(float))
    if still is not None:
        cov["dstill"] = pairwise_absdiff(still.astype(float))
    if turn is not None:
        cov["dturn"] = pairwise_absdiff(turn.astype(float))
    if pos_x is not None and pos_y is not None:
        cov["dspace_m"] = pairwise_euclid_xy(pos_x.astype(float), pos_y.astype(float))
    else:
        print('Position missing!')

    # NEW: distance-to-center difference (captures radial boundary context in circular arena)
    if dist_center is not None:
        cov["dcenter_m"] = pairwise_absdiff(dist_center.astype(float))

    # NEW: head-direction relative to center (egocentric bearing), circular distance
    if hd_rel_center_mean is not None:
        cov["dhd_rel_center_rad"] = circ_dist(
            hd_rel_center_mean.astype(float),
            hd_rel_center_mean.astype(float)
        )

    # ---- Validity mask for pairwise modeling
    valid = np.ones(n_ep, dtype=bool)
    if hd_valid is not None:
        valid &= (np.asarray(hd_valid, float) >= 0.7)  # threshold can be changed
    if pos_x is None or pos_y is None:
        # still can proceed; just no dspace
        pass

    tri_mask = upper_tri_mask_from_valid(valid)

    # ---- Vectorized forms for later partial correlations
    vec = {}
    vec["S_lat"] = upper_tri_vector(S_lat)
    if S_pv is not None:
        vec["S_pv"] = upper_tri_vector(S_pv)

    for k, M in cov.items():
        vec[k] = upper_tri_vector(M)

    vec["valid_pair_mask"] = tri_mask.astype(np.uint8)

    # ---- Save to HDF5
    meta = {
        "source_trajectories_h5": trajectories_h5,
        "source_behavior_h5": behavior_h5,
        "n_ep": int(n_ep),
        "win_bins": int(win_bins),
        "bin_size_s": float(bin_size_s),
        "strip_first_bins": int(strip_first_bins),
        "time_bins_used": [int(use_time_bins[0]), int(use_time_bins[-1]), int(use_time_bins.size)],
        "latent_similarity": latent_similarity,
        "pv_similarity": pv_similarity if S_pv is not None else None,
        "pv_transform": pv_transform if S_pv is not None else None,
        "pv_zscore": bool(pv_zscore) if S_pv is not None else None,
        "has_S_pv": bool(S_pv is not None),
        "has_dspace": bool("dspace_m" in cov),
        "valid_hd_thr": 0.7 if hd_valid is not None else None,
        "ep_select_idx_first_last": [int(ep_select_idx[0]), int(ep_select_idx[-1])] if len(ep_select_idx) else None,
        "note": "Stores episode similarity matrices and covariate matrices + upper-tri vectors for partial corr."
    }
    
    # sanity: new covariates should be finite for most pairs if present
    for k in ["dcenter_m", "dhd_rel_center_rad"]:
        if k in cov:
            frac = np.mean(np.isfinite(cov[k][np.triu_indices(n_ep, 1)]))
            print(f"[cov sanity] {k}: finite fraction={frac:.3f}")

    with h5py.File(out_h5_path, "a") as f:
        if group in f and overwrite:
            del f[group]
        g = f.require_group(group)
        g.attrs["meta_json"] = json.dumps(meta)

        # Similarities
        gs = g.require_group("similarity")
        gs.create_dataset("S_latent", data=S_lat.astype(np.float32), compression="gzip")
        if S_pv is not None:
            gs.create_dataset("S_pv", data=S_pv.astype(np.float32), compression="gzip")

        # Covariates
        gc = g.require_group("covariates")
        for k, M in cov.items():
            gc.create_dataset(k, data=M.astype(np.float32), compression="gzip")

        # Episode-level arrays (aligned)
        ge = g.require_group("episode_level")
        ge.create_dataset("valid_episode", data=valid.astype(np.uint8), compression="gzip")
        ge.create_dataset("target_bin_windows", data=target_bin_windows.astype(np.int64), compression="gzip")
        if hd_mean is not None:
            ge.create_dataset("hd_mean_rad", data=hd_mean.astype(np.float32), compression="gzip")
        if hd_valid is not None:
            ge.create_dataset("hd_valid_frac", data=hd_valid.astype(np.float32), compression="gzip")
        if still is not None:
            ge.create_dataset("stillness_mps", data=still.astype(np.float32), compression="gzip")
        if turn is not None:
            ge.create_dataset("head_turn_rate_rads", data=turn.astype(np.float32), compression="gzip")
        if pos_x is not None and pos_y is not None:
            ge.create_dataset("pos_mean_x_m", data=pos_x.astype(np.float32), compression="gzip")
            ge.create_dataset("pos_mean_y_m", data=pos_y.astype(np.float32), compression="gzip")

        # Vectorized upper-tri for modeling
        gv = g.require_group("vectors_uppertri")
        for k, v in vec.items():
            # scalars not expected here; all 1D
            gv.create_dataset(k, data=np.asarray(v), compression="gzip")

    return out_h5_path
