from __future__ import annotations

import os
import json
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import h5py
import numpy as np
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge


# ----------------------------
# Utilities
# ----------------------------

def select_window_starts_strided(
    valid_starts: np.ndarray,
    win_bins: int,
    stride_bins: int | None = None,
) -> np.ndarray:
    """
    Convert a dense list of all-valid sliding-window starts into a sparse list
    of starts, enforcing a minimum separation to avoid massive overlap.

    If stride_bins is None -> default stride = win_bins (non-overlapping).
    """
    if valid_starts.size == 0:
        return valid_starts

    if stride_bins is None:
        stride_bins = int(win_bins)

    valid_starts = np.asarray(valid_starts, dtype=np.int64)
    valid_starts.sort()

    chosen = []
    next_allowed = -10**18
    for s in valid_starts:
        if s >= next_allowed:
            chosen.append(int(s))
            next_allowed = int(s) + int(stride_bins)

    return np.asarray(chosen, dtype=np.int64)


def clip_periods(periods_bin: np.ndarray, t_bins: int) -> np.ndarray:
    """Clip periods to [0, t_bins-1] and drop invalid (end < start)."""
    p = np.asarray(periods_bin, dtype=np.int64).copy()
    if p.size == 0:
        return p.reshape(0, 2)
    p[:, 0] = np.clip(p[:, 0], 0, t_bins - 1)
    p[:, 1] = np.clip(p[:, 1], 0, t_bins - 1)
    keep = p[:, 1] >= p[:, 0]
    return p[keep]


def periods_to_mask(periods_bin: np.ndarray, t_bins: int) -> np.ndarray:
    """Create boolean mask length t_bins that is True for bins within any period."""
    mask = np.zeros(t_bins, dtype=bool)
    for s, e in np.asarray(periods_bin, dtype=np.int64):
        mask[s:e + 1] = True
    return mask


def pulse_periods_to_bin_periods(
    periods_pulse: np.ndarray,
    *,
    bins_per_pulse: int,
    end_inclusive: bool = True,
) -> np.ndarray:
    """
    Convert pulse-index periods (4 Hz) to bin-index periods.

    If pulse bins are [p*bpp, (p+1)*bpp - 1] for each pulse p,
    then period (p0, p1) maps to:
      start_bin = p0*bpp
      end_bin   = (p1+1)*bpp - 1  (if end_inclusive)  else p1*bpp - 1
    """
    P = np.asarray(periods_pulse, dtype=np.int64)
    if P.size == 0:
        return P.reshape(0, 2)
    if P.ndim != 2 or P.shape[1] != 2:
        raise ValueError("periods_pulse must be Nx2")
    p0 = P[:, 0]
    p1 = P[:, 1]
    start_bin = p0 * bins_per_pulse
    if end_inclusive:
        end_bin = (p1 + 1) * bins_per_pulse - 1
    else:
        end_bin = p1 * bins_per_pulse - 1
    return np.stack([start_bin, end_bin], axis=1).astype(np.int64)


def extract_fixed_windows_from_periods(
    periods_bin: np.ndarray,
    *,
    win_bins: int,
    t_bins: int,
    align: str = "start",
) -> np.ndarray:
    """
    From (start_bin, end_bin) target periods, extract fixed-length windows.
    align: "start" uses start_bin; "end" uses end_bin-win_bins+1; "center" centers window.
    Returns windows as Nx2 (start_bin, end_bin) inclusive.
    """
    P = np.asarray(periods_bin, dtype=np.int64)
    if P.size == 0:
        return P.reshape(0, 2)
    out = []
    for s, e in P:
        if align == "start":
            ws = s
        elif align == "end":
            ws = e - win_bins + 1
        elif align == "center":
            c = (s + e) // 2
            ws = c - win_bins // 2
        else:
            raise ValueError(f"Unknown align={align}")
        we = ws + win_bins - 1
        if ws < 0 or we >= t_bins:
            continue
        out.append((ws, we))
    return np.asarray(out, dtype=np.int64)


def preprocess_counts(
    X_counts: np.ndarray,
    *,
    X_is_neurons_by_time: bool,
    transform: str,
    fit_mask: np.ndarray,
) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
    """
    X_counts: either (n_neurons, t_bins) if X_is_neurons_by_time else (t_bins, n_neurons).
    Returns Xz (t_bins, n_neurons) z-scored using fit_mask; plus {'mu','sd'}.
    """
    X = X_counts.T if X_is_neurons_by_time else X_counts
    X = np.asarray(X, dtype=np.float64)
    if X.ndim != 2:
        raise ValueError("X_counts must be 2D")
    if transform == "sqrt":
        X = np.sqrt(np.clip(X, 0, None))
    elif transform == "log1p":
        X = np.log1p(np.clip(X, 0, None))
    elif transform == "none":
        pass
    else:
        raise ValueError(f"Unknown transform={transform}")

    fit_mask = np.asarray(fit_mask, dtype=bool)
    if fit_mask.shape[0] != X.shape[0]:
        raise ValueError("fit_mask length must match t_bins")

    mu = np.nanmean(X[fit_mask], axis=0)
    sd = np.nanstd(X[fit_mask], axis=0)
    sd = np.where(sd > 1e-8, sd, 1.0)
    Xz = (X - mu) / sd
    return Xz, {"mu": mu.astype(np.float32), "sd": sd.astype(np.float32)}


# ----------------------------
# Residualization helpers
# ----------------------------

def _bin_centers_from_edges(t_edges: Optional[np.ndarray], t_bins: int, bin_size_s: float) -> np.ndarray:
    if t_edges is not None:
        t_edges = np.asarray(t_edges, dtype=np.float64)
        if t_edges.ndim != 1 or t_edges.size != t_bins + 1:
            # fallback
            return (np.arange(t_bins, dtype=np.float64) + 0.5) * float(bin_size_s)
        return 0.5 * (t_edges[:-1] + t_edges[1:])
    return (np.arange(t_bins, dtype=np.float64) + 0.5) * float(bin_size_s)


def _interp_to_bin_centers(time_s_src: np.ndarray, y_src: np.ndarray, t_bin: np.ndarray) -> np.ndarray:
    """
    Linear interpolation for continuous covariates. Returns float64 vector len(t_bin).
    If y_src has NaNs, we interpolate on valid points only; if too few valid, return NaNs.
    """
    time_s_src = np.asarray(time_s_src, dtype=np.float64)
    y_src = np.asarray(y_src, dtype=np.float64)
    if y_src.ndim != 1:
        raise ValueError("y_src must be 1D")
    m = np.isfinite(time_s_src) & np.isfinite(y_src)
    if m.sum() < 2:
        return np.full_like(t_bin, np.nan, dtype=np.float64)
    return np.interp(t_bin, time_s_src[m], y_src[m], left=np.nan, right=np.nan)


def _unwrap_angle(angle_rad: np.ndarray) -> np.ndarray:
    a = np.asarray(angle_rad, dtype=np.float64)
    if a.size == 0:
        return a
    # unwrap on finite subset, keep NaNs
    out = a.copy()
    m = np.isfinite(out)
    if m.sum() >= 2:
        out[m] = np.unwrap(out[m])
    return out


def build_stim_phase_design_matrix(
    t_bins: int,
    *,
    bins_per_pulse: int,
    drop_first: bool = True,
    include_intercept: bool = True,
) -> Tuple[np.ndarray, List[str]]:
    """
    Build one-hot phase-of-pulse design:
      phase_idx = bin_index % bins_per_pulse
    With intercept, we drop_first to avoid collinearity.
    """
    if bins_per_pulse <= 1:
        raise ValueError("bins_per_pulse must be >= 2 for phase regression")

    phase = (np.arange(t_bins, dtype=np.int64) % int(bins_per_pulse)).astype(np.int64)
    n_levels = int(bins_per_pulse)

    # one-hot
    OH = np.zeros((t_bins, n_levels), dtype=np.float64)
    OH[np.arange(t_bins), phase] = 1.0

    names = [f"stim_phase_{k}" for k in range(n_levels)]
    if drop_first:
        OH = OH[:, 1:]
        names = names[1:]

    if include_intercept:
        X = np.column_stack([np.ones(t_bins, dtype=np.float64), OH])
        names = ["intercept"] + names
    else:
        X = OH
    return X, names


@dataclass
class ContextCovariateSpec:
    use_pos_xy: bool = True
    use_dist_center: bool = True
    use_hd_unwrap: bool = True
    use_head_ang_vel: bool = True
    use_body_speed: bool = True
    use_session_time: bool = True          # include t and t^2
    zscore_cols: bool = True               # zscore covariates (not intercept)


def build_context_design_matrix_for_bins(
    behavior_h5_path: str,
    *,
    t_edges_50ms: Optional[np.ndarray],
    t_bins: int,
    bin_size_s: float,
    cov_spec: ContextCovariateSpec = ContextCovariateSpec(),
) -> Tuple[np.ndarray, List[str]]:
    """
    Load dlc_features/session_ts covariates and interpolate them to neural bin centers.
    Requires group 'dlc_features/session_ts' to contain:
      - time_s
      - pos_x_m, pos_y_m
      - hd_angle_rad (wrapped) OR hd_angle_rad_unwrap (optional)
      - head_ang_vel_rads
      - body_speed_rms_mps
      - dist_center_m (optional)
    Returns Xcov with intercept in col0, and names.
    """
    with h5py.File(behavior_h5_path, "r") as f:
        if "dlc_features/session_ts" not in f:
            raise KeyError("behavior_h5 missing group dlc_features/session_ts")
        g = f["dlc_features/session_ts"]
        time_s = np.asarray(g["time_s"][...], dtype=np.float64)

        def req(name: str) -> np.ndarray:
            if name not in g:
                raise KeyError(f"behavior_h5 dlc_features/session_ts missing dataset {name}")
            return np.asarray(g[name][...], dtype=np.float64)

        pos_x = req("pos_x_m")
        pos_y = req("pos_y_m")

        if "hd_angle_rad_unwrap" in g:
            hd_u = np.asarray(g["hd_angle_rad_unwrap"][...], dtype=np.float64)
        else:
            hd = req("hd_angle_rad")
            hd_u = _unwrap_angle(hd)

        head_w = req("head_ang_vel_rads")
        sp = req("body_speed_rms_mps")

        dist_c = None
        if cov_spec.use_dist_center:
            if "dist_center_m" in g:
                dist_c = np.asarray(g["dist_center_m"][...], dtype=np.float64)
            else:
                # can compute from pos if missing
                dist_c = np.sqrt(pos_x**2 + pos_y**2)

    t_bin = _bin_centers_from_edges(t_edges_50ms, t_bins=t_bins, bin_size_s=bin_size_s)

    cols = []
    names = []

    if cov_spec.use_pos_xy:
        cols.append(_interp_to_bin_centers(time_s, pos_x, t_bin)); names.append("pos_x_m")
        cols.append(_interp_to_bin_centers(time_s, pos_y, t_bin)); names.append("pos_y_m")

    if cov_spec.use_dist_center and dist_c is not None:
        cols.append(_interp_to_bin_centers(time_s, dist_c, t_bin)); names.append("dist_center_m")

    if cov_spec.use_hd_unwrap:
        cols.append(_interp_to_bin_centers(time_s, hd_u, t_bin)); names.append("hd_angle_rad_unwrap")

    if cov_spec.use_head_ang_vel:
        cols.append(_interp_to_bin_centers(time_s, head_w, t_bin)); names.append("head_ang_vel_rads")

    if cov_spec.use_body_speed:
        cols.append(_interp_to_bin_centers(time_s, sp, t_bin)); names.append("body_speed_rms_mps")

    if cov_spec.use_session_time:
        # normalized time helps conditioning
        tt = (t_bin - np.nanmin(t_bin))
        tt = tt / (np.nanmax(tt) + 1e-9)
        cols.append(tt); names.append("session_t")
        cols.append(tt**2); names.append("session_t2")

    if len(cols) == 0:
        raise ValueError("No context covariates selected")

    X = np.column_stack(cols).astype(np.float64)

    # Build intercept + optionally z-score non-intercept columns
    m_all = np.isfinite(X).all(axis=1)
    if m_all.sum() < max(50, 0.1 * t_bins):
        # too many NaNs -> better to fail early than silently residualize garbage
        raise RuntimeError("Too few valid time bins for context design matrix (NaNs after interpolation).")

    if cov_spec.zscore_cols:
        mu = np.nanmean(X[m_all], axis=0)
        sd = np.nanstd(X[m_all], axis=0)
        sd = np.where(sd > 1e-8, sd, 1.0)
        X = (X - mu) / sd

    Xcov = np.column_stack([np.ones(t_bins, dtype=np.float64), X])
    names = ["intercept"] + names
    return Xcov, names


def residualize_latent_time_series(
    Z_all: np.ndarray,
    Xcov: np.ndarray,
    *,
    ridge_alpha: float = 1.0,
) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
    """
    Residualize each latent dimension independently using Ridge:
      z_d ~ Xcov
      z_res = z - z_hat
    Returns (Z_res, info) where info contains coef_ and intercept_ per dimension.
    """
    Z = np.asarray(Z_all, dtype=np.float64)
    X = np.asarray(Xcov, dtype=np.float64)
    if Z.ndim != 2 or X.ndim != 2:
        raise ValueError("Z_all and Xcov must be 2D")
    if Z.shape[0] != X.shape[0]:
        L = min(Z.shape[0], X.shape[0])
        Z = Z[:L]
        X = X[:L]

    m = np.isfinite(X).all(axis=1) & np.isfinite(Z).all(axis=1)
    if m.sum() < max(50, 0.1 * Z.shape[0]):
        raise RuntimeError("Too few valid rows for residualization")

    Zhat = np.zeros_like(Z, dtype=np.float64)
    coefs = np.zeros((Z.shape[1], X.shape[1]), dtype=np.float64)

    for d in range(Z.shape[1]):
        y = Z[m, d]
        Xd = X[m]
        model = Ridge(alpha=float(ridge_alpha), fit_intercept=False)
        model.fit(Xd, y)
        coefs[d, :] = model.coef_
        Zhat[m, d] = Xd @ model.coef_

    Zres = Z - Zhat
    info = {"coef": coefs.astype(np.float32)}
    return Zres.astype(np.float32), Zhat.astype(np.float32), info


def extract_episode_stack_from_Z(
    Z_all: np.ndarray,
    tgt_windows: np.ndarray,
    *,
    win_bins: int,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Stack episode windows from Z_all into traj_stack + traj_ptr.
    """
    Z_all = np.asarray(Z_all, dtype=np.float32)
    tgt_windows = np.asarray(tgt_windows, dtype=np.int64)
    n_ep = tgt_windows.shape[0]
    n_pc = Z_all.shape[1]

    traj_stack = np.empty((n_ep * win_bins, n_pc), dtype=np.float32)
    ep_ptr = np.zeros((n_ep, 2), dtype=np.int64)

    row = 0
    used = 0
    for i, (bs, be) in enumerate(tgt_windows):
        seg = Z_all[bs:be + 1]
        if seg.shape[0] != win_bins:
            continue
        traj_stack[row:row + win_bins] = seg
        ep_ptr[used] = (row, row + win_bins - 1)
        row += win_bins
        used += 1

    traj_stack = traj_stack[:used * win_bins]
    ep_ptr = ep_ptr[:used]
    return traj_stack, ep_ptr


def episode_mean_residualize_stack(
    traj_stack: np.ndarray,
    traj_ptr: np.ndarray,
):
    """
    Given a concatenated traj_stack (T_concat, D) and traj_ptr (n_ep, 2) inclusive indices into traj_stack,
    compute:
      - mu_ep: (n_ep, D) mean vector per episode (mean across time bins)
      - traj_stack_epmean: (T_concat, D) with episode mean subtracted within each episode

    Notes:
      - Does NOT change traj_ptr or episode membership.
      - Leaves NaNs in place, but uses nanmean so partial NaNs don't nuke everything.
    """
    if traj_stack.ndim != 2:
        raise ValueError(f"traj_stack must be 2D (T,D), got {traj_stack.shape}")
    if traj_ptr.ndim != 2 or traj_ptr.shape[1] != 2:
        raise ValueError(f"traj_ptr must be (n_ep,2), got {traj_ptr.shape}")

    n_ep = traj_ptr.shape[0]
    D = traj_stack.shape[1]

    mu_ep = np.full((n_ep, D), np.nan, dtype=np.float32)
    out = traj_stack.astype(np.float32, copy=True)

    for i in range(n_ep):
        s, e = int(traj_ptr[i, 0]), int(traj_ptr[i, 1])
        if s < 0 or e < s or e >= out.shape[0]:
            raise ValueError(f"Bad traj_ptr[{i}] = ({s},{e}) for stack len {out.shape[0]}")
        seg = out[s : e + 1, :]  # inclusive
        mu = np.nanmean(seg, axis=0)
        mu_ep[i] = mu
        out[s : e + 1, :] = seg - mu[None, :]

    return out, mu_ep


# ----------------------------
# Main builder
# ----------------------------

def build_core_trajectories_h5(
    out_h5_path: str,
    X_counts_50ms: np.ndarray,
    state_periods_pulse: dict,
    t_edges_50ms: np.ndarray | None = None,
    subgroup: str = None,
    mode: str = "w",
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
    align: str = "start",

    # residualization inputs
    behavior_h5_path: Optional[str] = None,
    do_context_residual: bool = True,
    do_stim_phase_residual: bool = True,
    context_cov_spec: ContextCovariateSpec = ContextCovariateSpec(),
    ridge_alpha_resid: float = 1.0,
):
    """
    Builds and stores:
      - stationary-fit z-scoring stats
      - PCA fit on stationary-only bins (bgr_sta + sil_sta, excluding target)
      - projections for all bins
      - extracted fixed-length target episode trajectories in PCA space

    PLUS (optional):
      - context-residualized latent time series + episodes (backward compatible as *_resid)
      - stim-phase residualized latent time series + episodes
      - context+stim residualized latent time series + episodes

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

    # Masks
    mask_tgt = periods_to_mask(tgt_bin, t_bins)
    mask_sta = periods_to_mask(bgr_bin, t_bins) | periods_to_mask(sil_bin, t_bins)
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
    if tgt_windows.shape[0] == 0:
        raise RuntimeError("No target windows survived (likely near edges). Adjust episode_win_s or alignment.")

    # Raw episode stack (baseline)
    traj_stack, ep_ptr = extract_episode_stack_from_Z(Z_all, tgt_windows, win_bins=win_bins)

    # IMPORTANT: ep_binwin should match used episodes (in case any were skipped)
    used_ep = ep_ptr.shape[0]
    ep_binwin = tgt_windows[:used_ep].copy()

    # Residualization outputs (optional)
    Z_all_resid = None
    traj_stack_resid = None
    Z_all_resid_stim = None
    traj_stack_resid_stim = None
    Z_all_resid_ctx_stim = None
    traj_stack_resid_ctx_stim = None

    resid_meta = {}

    # Build stim-phase design once (always possible)
    Xstim, stim_names = build_stim_phase_design_matrix(
        t_bins=t_bins,
        bins_per_pulse=bins_per_pulse,
        drop_first=True,
        include_intercept=True,
    )

    # Context design (requires behavior file)
    Xctx = None
    ctx_names = None
    if behavior_h5_path is not None and do_context_residual:
        #context_cov_spec.use_session_time = False  # disable session time for residualization
        #context_cov_spec.use_pos_xy = False
        Xctx, ctx_names = build_context_design_matrix_for_bins(
            behavior_h5_path,
            t_edges_50ms=t_edges_50ms,
            t_bins=t_bins,
            bin_size_s=bin_size_s,
            cov_spec=context_cov_spec,
        )

    # (A) Context residual (backward compatible: *_resid)
    if Xctx is not None:
        Z_all_resid, Zhat_resid, info_ctx = residualize_latent_time_series(Z_all, Xctx, ridge_alpha=ridge_alpha_resid)
        traj_stack_resid, _ = extract_episode_stack_from_Z(Z_all_resid, ep_binwin, win_bins=win_bins)
        traj_stack_hat_resid, _ = extract_episode_stack_from_Z(Zhat_resid, ep_binwin, win_bins=win_bins)
        resid_meta["context"] = {
            "ridge_alpha": float(ridge_alpha_resid),
            "covariates": list(ctx_names),
            "coef_shape": list(info_ctx["coef"].shape),
        }

    # (B) Stim-phase residual
    if do_stim_phase_residual:
        Z_all_resid_stim, Zhat_resid_stim, info_stim = residualize_latent_time_series(Z_all, Xstim, ridge_alpha=ridge_alpha_resid)
        traj_stack_resid_stim, _ = extract_episode_stack_from_Z(Z_all_resid_stim, ep_binwin, win_bins=win_bins)
        traj_stack_hat_resid_stim, _ = extract_episode_stack_from_Z(Zhat_resid_stim, ep_binwin, win_bins=win_bins)
        resid_meta["stim_phase"] = {
            "ridge_alpha": float(ridge_alpha_resid),
            "covariates": list(stim_names),
            "coef_shape": list(info_stim["coef"].shape),
        }

    # (C) Context + Stim combined residual
    if (Xctx is not None) and do_stim_phase_residual:
        # both have intercept in col0 -> remove one intercept to avoid duplication
        Xboth = np.column_stack([Xctx, Xstim[:, 1:]])  # keep ctx intercept, drop stim intercept
        both_names = list(ctx_names) + [n for n in stim_names if n != "intercept"]
        Z_all_resid_ctx_stim, Zhat_resid_ctx_stim, info_both = residualize_latent_time_series(Z_all, Xboth, ridge_alpha=ridge_alpha_resid)
        traj_stack_resid_ctx_stim, _ = extract_episode_stack_from_Z(Z_all_resid_ctx_stim, ep_binwin, win_bins=win_bins)
        traj_stack_hat_resid_ctx_stim, _ = extract_episode_stack_from_Z(Zhat_resid_ctx_stim, ep_binwin, win_bins=win_bins)
        resid_meta["context_plus_stim"] = {
            "ridge_alpha": float(ridge_alpha_resid),
            "covariates": list(both_names),
            "coef_shape": list(info_both["coef"].shape),
        }

    # ---------------- Episode-mean residualization (within-episode mean subtraction) ----------------
    # Raw
    traj_stack_epmean, traj_mu = episode_mean_residualize_stack(traj_stack, ep_ptr)

    # Context residualized (back-compat name: traj_stack_resid)
    traj_stack_resid_epmean, traj_mu_resid = episode_mean_residualize_stack(traj_stack_resid, ep_ptr)
    traj_stack_hat_resid_epmean, traj_mu_hat_resid = episode_mean_residualize_stack(traj_stack_hat_resid, ep_ptr)

    # Stim-only residualized
    traj_stack_resid_stim_epmean, traj_mu_resid_stim = episode_mean_residualize_stack(traj_stack_resid_stim, ep_ptr)
    traj_stack_hat_resid_stim_epmean, traj_mu_hat_resid_stim = episode_mean_residualize_stack(traj_stack_hat_resid_stim, ep_ptr)

    # Context+stim residualized
    traj_stack_resid_ctx_stim_epmean, traj_mu_resid_ctx_stim = episode_mean_residualize_stack(traj_stack_resid_ctx_stim, ep_ptr)
    traj_stack_hat_resid_ctx_stim_epmean, traj_mu_hat_resid_ctx_stim = episode_mean_residualize_stack(traj_stack_hat_resid_ctx_stim, ep_ptr)

    # Meta
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
        residualization=dict(
            behavior_h5_path=str(behavior_h5_path) if behavior_h5_path is not None else "",
            do_context_residual=bool(do_context_residual),
            do_stim_phase_residual=bool(do_stim_phase_residual),
            ridge_alpha=float(ridge_alpha_resid),
        ),
    )
    meta["episode_mean_residualization"] = True
    meta["episode_mean_residualization_outputs"] = {
        "traj_stack_epmean": "episodes/traj_stack_epmean",
        "traj_mu": "episodes/traj_mu",
        "traj_stack_resid_epmean": "episodes/traj_stack_resid_epmean",
        "traj_mu_resid": "episodes/traj_mu_resid",
        "traj_stack_resid_stim_epmean": "episodes/traj_stack_resid_stim_epmean",
        "traj_mu_resid_stim": "episodes/traj_mu_resid_stim",
        "traj_stack_resid_ctx_stim_epmean": "episodes/traj_stack_resid_ctx_stim_epmean",
        "traj_mu_resid_ctx_stim": "episodes/traj_mu_resid_ctx_stim",
    }

    # Write HDF5
    with h5py.File(out_h5_path, mode) as hfile:
        if subgroup is not None:
            f = hfile.create_group(subgroup)
        else:
            f = hfile

        f.attrs["meta_json"] = json.dumps(meta)

        g_in = f.create_group("inputs")
        if t_edges_50ms is not None:
            g_in.create_dataset("t_edges_50ms", data=np.asarray(t_edges_50ms, dtype=np.float64), compression="gzip")

        # states (bin units)
        g_state = f.create_group("states")
        g_state.create_dataset("target_periods_bin", data=tgt_bin, compression="gzip")
        g_state.create_dataset("bgr_sta_periods_bin", data=bgr_bin, compression="gzip")
        g_state.create_dataset("sil_sta_periods_bin", data=sil_bin, compression="gzip")
        g_state.create_dataset("mask_tgt", data=mask_tgt.astype(np.uint8), compression="gzip")
        g_state.create_dataset("mask_sta", data=mask_sta.astype(np.uint8), compression="gzip")
        g_state.create_dataset("mask_pca_fit", data=mask_pca_fit.astype(np.uint8), compression="gzip")

        # preprocessing stats
        g_pp = f.create_group("preproc")
        g_pp.create_dataset("z_mu", data=zstats["mu"], compression="gzip")
        g_pp.create_dataset("z_sd", data=zstats["sd"], compression="gzip")

        # PCA model
        g_pca = f.create_group("pca")
        g_pca.create_dataset("components", data=pca.components_.astype(np.float32), compression="gzip")
        g_pca.create_dataset("mean", data=pca.mean_.astype(np.float32), compression="gzip")
        g_pca.create_dataset("explained_variance", data=pca.explained_variance_.astype(np.float32), compression="gzip")
        g_pca.create_dataset("explained_variance_ratio", data=pca.explained_variance_ratio_.astype(np.float32), compression="gzip")

        # latent projections
        g_lat = f.create_group("latent")
        g_lat.create_dataset("Z_all", data=Z_all, compression="gzip")

        if Z_all_resid is not None:
            g_lat.create_dataset("Z_all_resid", data=Z_all_resid, compression="gzip")
            g_lat.create_dataset("Zhat_resid", data=Zhat_resid, compression="gzip")
        if Z_all_resid_stim is not None:
            g_lat.create_dataset("Z_all_resid_stim", data=Z_all_resid_stim, compression="gzip")
            g_lat.create_dataset("Zhat_resid_stim", data=Zhat_resid_stim, compression="gzip")
        if Z_all_resid_ctx_stim is not None:
            g_lat.create_dataset("Z_all_resid_ctx_stim", data=Z_all_resid_ctx_stim, compression="gzip")
            g_lat.create_dataset("Zhat_resid_ctx_stim", data=Zhat_resid_ctx_stim, compression="gzip")

        # episodes
        g_ep = f.create_group("episodes")
        g_ep.create_dataset("target_bin_windows", data=ep_binwin, compression="gzip")
        g_ep.create_dataset("traj_ptr", data=ep_ptr, compression="gzip")
        g_ep.create_dataset("traj_stack", data=traj_stack, compression="gzip")

        if traj_stack_resid is not None:
            g_ep.create_dataset("traj_stack_resid", data=traj_stack_resid, compression="gzip")
            g_ep.create_dataset("traj_stack_hat_resid", data=traj_stack_hat_resid, compression="gzip")
        if traj_stack_resid_stim is not None:
            g_ep.create_dataset("traj_stack_resid_stim", data=traj_stack_resid_stim, compression="gzip")
            g_ep.create_dataset("traj_stack_hat_resid_stim", data=traj_stack_hat_resid_stim, compression="gzip")
        if traj_stack_resid_ctx_stim is not None:
            g_ep.create_dataset("traj_stack_resid_ctx_stim", data=traj_stack_resid_ctx_stim, compression="gzip")
            g_ep.create_dataset("traj_stack_hat_resid_ctx_stim", data=traj_stack_hat_resid_ctx_stim, compression="gzip")

        mu_dict = {
            "traj_stack_epmean": traj_stack_epmean,
            "traj_mu": traj_mu,

            "traj_stack_resid_epmean": traj_stack_resid_epmean,
            "traj_mu_resid": traj_mu_resid,
            "traj_stack_resid_epmean": traj_stack_hat_resid_epmean,
            "traj_mu_resid": traj_mu_hat_resid,

            "traj_stack_resid_stim_epmean": traj_stack_resid_stim_epmean,
            "traj_mu_resid_stim": traj_mu_resid_stim,
            "traj_stack_resid_stim_epmean": traj_stack_hat_resid_stim_epmean,
            "traj_mu_resid_stim": traj_mu_hat_resid_stim,

            "traj_stack_resid_ctx_stim_epmean": traj_stack_resid_ctx_stim_epmean,
            "traj_mu_resid_ctx_stim": traj_mu_resid_ctx_stim,
            "traj_stack_resid_ctx_stim_epmean": traj_stack_hat_resid_ctx_stim_epmean,
            "traj_mu_resid_ctx_stim": traj_mu_hat_resid_ctx_stim,

        }
        for k, v in mu_dict.items():
            # if k in g_ep:
            #     del g_ep[k]
            g_ep.create_dataset(k, data=v, compression="gzip")

        # residualization metadata
        g_r = f.create_group("residualization")
        g_r.attrs["meta_json"] = json.dumps(resid_meta)

        # store stim design info (for debugging)
        g_r.create_dataset("stim_phase_names", data=np.array(stim_names, dtype="S"), compression="gzip")
        if ctx_names is not None:
            g_r.create_dataset("context_names", data=np.array(ctx_names, dtype="S"), compression="gzip")

    return out_h5_path
