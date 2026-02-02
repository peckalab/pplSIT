import json
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
import h5py

from sklearn.linear_model import PoissonRegressor, Ridge
from sklearn.model_selection import KFold


# -----------------------------
# Utilities
# -----------------------------

def _windows_to_mask(windows: np.ndarray, T: int) -> np.ndarray:
    """
    Convert (n_win, 2) [start,end] (inclusive) windows to boolean mask of length T.
    """
    m = np.zeros(T, dtype=bool)
    if windows.size == 0:
        return m
    w = np.asarray(windows, dtype=np.int64)
    if w.ndim != 2 or w.shape[1] != 2:
        raise ValueError("windows must have shape (n,2)")
    for s, e in w:
        s = int(max(0, s))
        e = int(min(T - 1, e))
        if e >= s:
            m[s : e + 1] = True
    return m


def _load_fit_mask_from_trajectories(
    trajectories_h5_path: str,
    T: int,
    fit_period: str,
) -> Optional[np.ndarray]:
    """
    fit_period:
      - "all":     no masking (return None)
      - "target":  states/mask_tgt
      - "all_sta": states/mask_sta  (or windows in episodes_sta/all/bin_windows)
      - "bgr_sta": states/mask_bgr_sta OR states/bgr_sta_periods_bin OR episodes_sta/bgr/bin_windows
      - "sil_sta": states/mask_sil_sta OR states/sil_sta_periods_bin OR episodes_sta/sil/bin_windows
    Returns:
      mask (T,) bool OR None for "all".
    """
    fit_period = str(fit_period).lower()
    if fit_period in ("all", "full", "session"):
        return None

    with h5py.File(trajectories_h5_path, "r") as f:
        # all sta is mask_sta | mask_tgt
        if fit_period == "all_sta":
            m = f["states/mask_tgt"][...].astype(bool) | f["states/mask_sta"][...].astype(bool)
            if m.shape[0] != T:
                raise ValueError(f"mask length {m.shape[0]} != T {T}")
            return m

        # --- direct boolean masks (preferred)
        direct_mask_paths = {
            "target": ["states/mask_tgt"],
            "all_sta": ["states/mask_sta"],
            "bgr_sta": ["states/mask_bgr_sta"],
            "sil_sta": ["states/mask_sil_sta"],
        }
        for p in direct_mask_paths.get(fit_period, []):
            if p in f:
                m = f[p][...].astype(bool)
                if m.shape[0] != T:
                    raise ValueError(f"{p} length {m.shape[0]} != T {T}")
                return m

        # --- alternative: windows stored as (n,2)
        window_paths = []
        if fit_period == "all_sta":
            window_paths = [
                "episodes_sta/all/raw/bin_windows",
                "episodes_sta/all/bin_windows",
            ]
        elif fit_period == "bgr_sta":
            window_paths = [
                "states/bgr_sta_periods_bin",
                "episodes_sta/bgr/raw/bin_windows",
                "episodes_sta/bgr/bin_windows",
            ]
        elif fit_period == "sil_sta":
            window_paths = [
                "states/sil_sta_periods_bin",
                "episodes_sta/sil/raw/bin_windows",
                "episodes_sta/sil/bin_windows",
            ]
        elif fit_period == "target":
            window_paths = [
                "episodes/target_bin_windows",
                "episodes/bin_windows",
            ]

        for p in window_paths:
            if p in f:
                w = f[p][...]
                # if it's already a mask, accept it
                if w.ndim == 1 and w.shape[0] == T:
                    return w.astype(bool)
                return _windows_to_mask(w, T)

    raise KeyError(
        f"Could not find mask/windows for fit_period='{fit_period}' in {trajectories_h5_path}"
    )


def _safe_zscore(x: np.ndarray, eps: float = 1e-8) -> np.ndarray:
    x = np.asarray(x, dtype=np.float32)
    m = np.nanmean(x)
    s = np.nanstd(x)
    if not np.isfinite(s) or s < eps:
        return np.zeros_like(x, dtype=np.float32)
    return (x - m) / (s + eps)


def _interp_to_bins(t_src: np.ndarray, x_src: np.ndarray, t_bin: np.ndarray) -> np.ndarray:
    """
    Linear interpolate x_src defined at times t_src to bin-center times t_bin.
    x_src can be 1D or 2D (T_src, K). Returns (T_bin, K).
    """
    t_src = np.asarray(t_src).astype(np.float64)
    t_bin = np.asarray(t_bin).astype(np.float64)

    if x_src.ndim == 1:
        x_src = x_src[:, None]
    x_src = np.asarray(x_src, dtype=np.float64)

    K = x_src.shape[1]
    out = np.empty((len(t_bin), K), dtype=np.float64)

    # np.interp works 1D; do per column
    for k in range(K):
        y = x_src[:, k]
        # fill NaNs by nearest valid before interpolate
        ok = np.isfinite(y) & np.isfinite(t_src)
        if ok.sum() < 2:
            out[:, k] = 0.0
            continue
        out[:, k] = np.interp(t_bin, t_src[ok], y[ok], left=y[ok][0], right=y[ok][-1])
    return out.astype(np.float32)


def _bin_centers_from_edges(bin_edges_s: np.ndarray) -> np.ndarray:
    bin_edges_s = np.asarray(bin_edges_s, dtype=np.float64)
    return 0.5 * (bin_edges_s[:-1] + bin_edges_s[1:])


def _make_target_mask_and_phase(
    T: int,
    bin_centers_s: np.ndarray,
    target_bin_windows: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Build:
      mask_tgt: (T,) bool
      phi:      (T,) float in [0,1], only meaningful where mask_tgt True, else 0
    target_bin_windows is (n_ep, 2) in BIN INDICES [start, end_exclusive) or [start,end_inclusive]?
    We handle both by normalizing to half-open [s, e).
    """
    mask = np.zeros(T, dtype=bool)
    phi = np.zeros(T, dtype=np.float32)

    if target_bin_windows is None or len(target_bin_windows) == 0:
        return mask, phi

    w = np.asarray(target_bin_windows, dtype=np.int64)
    if w.ndim != 2 or w.shape[1] != 2:
        raise ValueError("target_bin_windows must be (n,2)")

    # detect inclusive vs exclusive: if any end==start, or ends exceed T, etc.
    # We'll convert each window to half-open [s, e) safely:
    for s, e in w:
        s = int(s)
        e = int(e)
        if e <= s:
            continue

        # If you stored inclusive end, then typical length would be e-s+1.
        # We can't know for sure; heuristic: if e == T-1 sometimes, it's probably inclusive.
        # We'll treat as:
        #   if e < T and (e - s) >= 1 and (e - s) <= (T-1): assume inclusive? -> make e = e+1
        # but if already exclusive, e might equal T. handle both.
        if e <= T - 1:
            # could be inclusive
            e_half_open = e + 1
        else:
            # could already be exclusive
            e_half_open = e

        e_half_open = max(s + 1, min(T, e_half_open))
        s = max(0, min(T - 1, s))

        idx = np.arange(s, e_half_open, dtype=np.int64)
        if len(idx) == 0:
            continue

        mask[idx] = True
        # phi from 0 to 1 across the window
        L = len(idx)
        phi[idx] = np.linspace(0.0, 1.0, L, endpoint=True).astype(np.float32)

    return mask, phi


def _make_stim_regressors(
    bin_centers_s: np.ndarray,
    stim_times_s: np.ndarray,
    f_hz: float = 4.0,
    impulse_kernel_s: float = 0.5,
    bin_size_s: float = 0.05,
) -> Tuple[np.ndarray, List[str]]:
    """
    Stim block:
      - sin/cos at f_hz using absolute time
      - impulse regressor: stimulus onsets binned + convolved with exp kernel
    """
    t = np.asarray(bin_centers_s, dtype=np.float64)
    T = len(t)

    # sin/cos
    ang = 2.0 * np.pi * f_hz * t
    sin_f = np.sin(ang).astype(np.float32)
    cos_f = np.cos(ang).astype(np.float32)

    # impulse train
    stim_times_s = np.asarray(stim_times_s, dtype=np.float64)
    imp = np.zeros(T, dtype=np.float32)
    if stim_times_s.size > 0:
        # map each stim time to nearest bin
        idx = np.searchsorted(t, stim_times_s, side="left")
        idx = np.clip(idx, 0, T - 1)
        imp[idx] += 1.0

    # convolve with exponential kernel (causal)
    tau = float(impulse_kernel_s)
    if tau > 0:
        L = int(np.ceil(5 * tau / bin_size_s))
        tt = np.arange(L, dtype=np.float64) * bin_size_s
        k = np.exp(-tt / tau)
        k /= (k.sum() + 1e-12)
        imp_f = np.convolve(imp, k.astype(np.float32), mode="full")[:T].astype(np.float32)
    else:
        imp_f = imp

    X = np.stack([sin_f, cos_f, imp_f], axis=1)
    names = [f"stim_sin_{f_hz:.3g}Hz", f"stim_cos_{f_hz:.3g}Hz", f"stim_imp_exp{impulse_kernel_s:.3g}s"]
    return X, names


def _pseudo_r2_poisson(y: np.ndarray, mu: np.ndarray, eps: float = 1e-12) -> float:
    """
    McFadden-ish pseudo-R2 using deviance ratio:
      R2 = 1 - Dev(model)/Dev(null)
    with Poisson deviance:
      Dev = 2 * sum( y*log(y/mu) - (y - mu) ), with convention y*log(y/mu)=0 if y==0.
    """
    y = np.asarray(y, dtype=np.float64)
    mu = np.asarray(mu, dtype=np.float64)
    mu = np.clip(mu, eps, None)

    ybar = np.mean(y)
    mu0 = np.full_like(y, max(ybar, eps), dtype=np.float64)

    def dev(yv, muv):
        term = np.zeros_like(yv)
        nz = yv > 0
        term[nz] = yv[nz] * np.log(yv[nz] / muv[nz])
        return 2.0 * np.sum(term - (yv - muv))

    dev_m = dev(y, mu)
    dev_0 = dev(y, mu0)
    if dev_0 <= eps:
        return 0.0
    return float(1.0 - dev_m / dev_0)


# -----------------------------
# Design specs
# -----------------------------

@dataclass
class GLMSpec:
    # model choice
    model: str = "poisson"  # "poisson" or "ridge"
    alpha: float = 1.0      # L2 strength for PoissonRegressor or Ridge

    # CV
    n_splits: int = 5
    random_seed: int = 0

    # Blocks: parameters
    drift_pcs_k: int = 3

    behavior_keys: Tuple[str, ...] = (
        "pos_x_m", "pos_y_m",
        "dist_to_center_m",
        "hd_angle_rad",
        "head_ang_vel_rads",
        "body_speed_rms_mps",
    )

    stim_f_hz: float = 4.0
    stim_impulse_kernel_s: float = 0.5

    include_phi2: bool = True


# -----------------------------
# Loaders for external HDF5 files
# -----------------------------

def load_behavior_to_bins(
    behavior_h5_path: str,
    bin_centers_s: np.ndarray,
    keys: Tuple[str, ...],
    session_ts_group: str = "session_ts",
) -> Tuple[np.ndarray, List[str]]:
    """
    behavior_features.h5 is expected to contain:
      /session_ts/time_s (seconds) and datasets for each key (same length)
    If your structure differs, adjust here.
    """
    with h5py.File(behavior_h5_path, "r") as hfile:
        f = hfile['dlc_features']
        if session_ts_group not in f:
            raise KeyError(f"Expected group '{session_ts_group}' in {behavior_h5_path}")
        g = f[session_ts_group]

        # timebase
        if "time_s" not in g:
            raise KeyError(f"Expected '{session_ts_group}/time_s' in {behavior_h5_path}")
        t_src = g["time_s"][...].astype(np.float64)

        X_cols = []
        names = []
        for k in keys:
            if k not in g:
                raise KeyError(f"Missing behavior key '{session_ts_group}/{k}' in {behavior_h5_path}")
            v = g[k][...].astype(np.float32)
            v_i = _interp_to_bins(t_src, v, bin_centers_s)[:, 0]
            v_i = _safe_zscore(v_i)
            X_cols.append(v_i)
            names.append(k)

    X = np.stack(X_cols, axis=1).astype(np.float32) if len(X_cols) else np.zeros((len(bin_centers_s), 0), np.float32)
    return X, names


def load_drift_pcs_to_bins(
    drift_h5_path: str,
    bin_centers_s: np.ndarray,
    condition: str = "target",
    k: int = 3,
) -> Tuple[np.ndarray, List[str]]:
    """
    session_drift_pcs.h5 expected structure:
      /<condition>/pca/scores  (n_ep, n_pc)
      /<condition>/episode_center_s  (n_ep,) episode time (seconds)
    We interpolate pc_scores over episode time to all bins.

    condition in {"target","bgr_sta","sil_sta","all_sta"} (depending on what you stored).
    """
    with h5py.File(drift_h5_path, "r") as hfile:
        f = hfile['drift_pcs']
        if condition not in f:
            raise KeyError(f"Condition '{condition}' not found in {drift_h5_path}")
        g = f[condition]
        if "scores" not in g['pca'] or "episode_center_s" not in g:
            raise KeyError(f"Expected datasets '{condition}/pca/scores' and '{condition}/episode_center_s' in {drift_h5_path}")

        scores = g["pca"]["scores"][...].astype(np.float32)
        ep_t = g["episode_center_s"][...].astype(np.float64)

    k = int(min(k, scores.shape[1]))
    if k <= 0:
        return np.zeros((len(bin_centers_s), 0), np.float32), []

    # interpolate each PC score over time -> bins
    X = _interp_to_bins(ep_t, scores[:, :k], bin_centers_s)
    # z-score each column for numerical stability
    for j in range(X.shape[1]):
        X[:, j] = _safe_zscore(X[:, j])

    names = [f"drift_pc{j+1}" for j in range(k)]
    return X.astype(np.float32), names


def load_target_windows_from_trajectories(
    trajectories_h5_path: str
) -> Tuple[np.ndarray, Dict]:
    """
    trajectories.h5 expected to store target windows in bins.
    We'll look in a few common places. Adjust if needed.

    Returns:
      target_bin_windows: (n_ep,2)
      meta: dict with bin_size_s if present
    """
    with h5py.File(trajectories_h5_path, "r") as f:
        meta = json.loads(f.attrs.get("meta_json", "{}")) if "meta_json" in f.attrs else {}
        # Prefer explicit
        cand_paths = [
            "episodes/target_bin_windows",
            "episodes/bin_windows",
            "episodes/target_windows_bins",
        ]
        w = None
        for p in cand_paths:
            if p in f:
                w = f[p][...].astype(np.int64)
                break
        if w is None:
            raise KeyError(f"Could not find target windows in {trajectories_h5_path} (tried {cand_paths})")
    return w, meta


# -----------------------------
# Fitters
# -----------------------------

def fit_unit_poisson_cv(
    X: np.ndarray,
    y: np.ndarray,
    alpha: float,
    n_splits: int,
    seed: int,
) -> Tuple[float, np.ndarray, float]:
    """
    Returns:
      r2_cv (mean across folds),
      coef (P,),
      intercept
    Coef/intercept from fit on full data (not fold-average).
    """
    y = np.asarray(y, dtype=np.float64)
    X = np.asarray(X, dtype=np.float64)

    kf = KFold(n_splits=n_splits, shuffle=True, random_state=seed)
    r2s = []
    for tr, te in kf.split(X):
        mdl = PoissonRegressor(alpha=alpha, max_iter=1000)
        mdl.fit(X[tr], y[tr])
        mu = mdl.predict(X[te])
        r2s.append(_pseudo_r2_poisson(y[te], mu))
    r2_cv = float(np.mean(r2s)) if len(r2s) else 0.0

    mdl = PoissonRegressor(alpha=alpha, max_iter=1000)
    mdl.fit(X, y)
    coef = mdl.coef_.astype(np.float32)
    intercept = float(mdl.intercept_)
    return r2_cv, coef, intercept


def fit_unit_ridge_cv(
    X: np.ndarray,
    y: np.ndarray,
    alpha: float,
    n_splits: int,
    seed: int,
) -> Tuple[float, np.ndarray, float]:
    """
    Ridge on z-scored y (or raw y; you choose upstream).
    Returns CV R2 (standard), coef/intercept from full fit.
    """
    y = np.asarray(y, dtype=np.float64)
    X = np.asarray(X, dtype=np.float64)

    kf = KFold(n_splits=n_splits, shuffle=True, random_state=seed)
    r2s = []
    for tr, te in kf.split(X):
        mdl = Ridge(alpha=alpha)
        mdl.fit(X[tr], y[tr])
        r2s.append(float(mdl.score(X[te], y[te])))
    r2_cv = float(np.mean(r2s)) if len(r2s) else 0.0

    mdl = Ridge(alpha=alpha)
    mdl.fit(X, y)
    coef = mdl.coef_.astype(np.float32)
    intercept = float(mdl.intercept_)
    return r2_cv, coef, intercept


# -----------------------------
# Main entry point
# -----------------------------

def compute_and_save_cell_glm_single_session(
    out_h5_path: str,
    Y_counts: np.ndarray,                 # (T, N)
    bin_edges_s: np.ndarray,              # (T+1,)
    bin_size_s: float,
    behavior_h5_path: str,
    stim_times_s: np.ndarray,
    trajectories_h5_path: str,
    drift_pcs_h5_path: str,
    drift_condition: str = "target",      # which drift series to use as "slow offset" regressor timebase
    fit_period: str = "all",
    spec: Optional[GLMSpec] = None,
    unit_meta: Optional[np.ndarray] = None,  # (N, M) optional placeholder
) -> None:
    """
    Builds X=[Intercept | A drift PCs | B behavior | C stim | D target phase] and fits per-unit GLM.
    Stores:
      - design names + block slices
      - full-model coefs/intercepts
      - CV pseudo-R2 per unit
      - drop-block delta-R2 (fit reduced model without that block)
    """
    if spec is None:
        spec = GLMSpec()

    Y = np.asarray(Y_counts, dtype=np.float32)
    if Y.ndim != 2:
        raise ValueError("Y_counts must be (T_bins, N_units)")
    T, N = Y.shape

    bin_centers_s = _bin_centers_from_edges(bin_edges_s)
    if len(bin_centers_s) != T:
        raise ValueError(f"bin_edges_s length {len(bin_edges_s)} inconsistent with Y T={T}")

    fit_mask = _load_fit_mask_from_trajectories(
        trajectories_h5_path=trajectories_h5_path,
        T=T,
        fit_period=fit_period,
    )
    print(f'MASKING: {fit_mask.sum()}' if fit_mask is not None else 'No masking applied!!!')

    # ---- Block A: drift PCs
    XA, namesA = load_drift_pcs_to_bins(
        drift_h5_path=drift_pcs_h5_path,
        bin_centers_s=bin_centers_s,
        condition=drift_condition,
        k=spec.drift_pcs_k,
    )

    # ---- Block B: behavior
    XB, namesB = load_behavior_to_bins(
        behavior_h5_path=behavior_h5_path,
        bin_centers_s=bin_centers_s,
        keys=spec.behavior_keys,
        session_ts_group="session_ts",
    )

    # ---- Block C: stim
    XC, namesC = _make_stim_regressors(
        bin_centers_s=bin_centers_s,
        stim_times_s=np.asarray(stim_times_s, dtype=np.float64),
        f_hz=spec.stim_f_hz,
        impulse_kernel_s=spec.stim_impulse_kernel_s,
        bin_size_s=bin_size_s,
    )
    # z-score stim columns too
    for j in range(XC.shape[1]):
        XC[:, j] = _safe_zscore(XC[:, j])

    # ---- Block D: target phase/progress
    target_windows, meta_traj = load_target_windows_from_trajectories(trajectories_h5_path)
    mask_tgt, phi = _make_target_mask_and_phase(T, bin_centers_s, target_windows)

    colsD = [ _safe_zscore((phi * mask_tgt.astype(np.float32))) ]
    namesD = ["tgt_phi"]
    if spec.include_phi2:
        colsD.append(_safe_zscore(((phi ** 2) * mask_tgt.astype(np.float32))))
        namesD.append("tgt_phi2")
    XD = np.stack(colsD, axis=1).astype(np.float32)

    # ---- Assemble full design
    # Intercept will be handled by sklearn intercept; we store names without intercept
    X_blocks = {"A_drift": XA, "B_behavior": XB, "C_stim": XC, "D_phase": XD}
    names_blocks = {"A_drift": namesA, "B_behavior": namesB, "C_stim": namesC, "D_phase": namesD}

    X = np.concatenate([XA, XB, XC, XD], axis=1).astype(np.float32)
    names = namesA + namesB + namesC + namesD
    P = X.shape[1]

    # record block slices in the concatenated X
    slices = {}
    c0 = 0
    for bn in ["A_drift", "B_behavior", "C_stim", "D_phase"]:
        k = X_blocks[bn].shape[1]
        slices[bn] = (c0, c0 + k)
        c0 += k

    # ---- Optionally restrict to finite rows
    good_row = np.isfinite(X).all(axis=1)
    # Also drop rows where all units are NaN (shouldn't happen) but keep it safe:
    good_row &= np.isfinite(Y).all(axis=1)

    Xg = X[good_row]
    Yg = Y[good_row]
    Tg = Xg.shape[0]
    if Tg < 100:
        raise ValueError(f"Too few valid time bins after cleaning: {Tg}")

    # ---- Fit per unit
    r2_full = np.zeros(N, dtype=np.float32)
    coef_full = np.zeros((N, P), dtype=np.float32)
    intercept_full = np.zeros(N, dtype=np.float32)

    # reduced fits for block-drop deltas
    delta_r2 = {bn: np.zeros(N, dtype=np.float32) for bn in slices.keys()}

    # Precompute reduced X matrices for efficiency
    reduced_X = {}
    for bn, (a, b) in slices.items():
        keep = np.r_[0:a, b:P]
        reduced_X[bn] = Xg[:, keep]

    # choose fitter
    # def fit_one(Xmat, yvec):
    #     if spec.model == "poisson":
    #         return fit_unit_poisson_cv(Xmat, yvec, alpha=spec.alpha, n_splits=spec.n_splits, seed=spec.random_seed)
    #     elif spec.model == "ridge":
    #         # common choice: z-score y for ridge
    #         return fit_unit_ridge_cv(Xmat, _safe_zscore(yvec), alpha=spec.alpha, n_splits=spec.n_splits, seed=spec.random_seed)
    #     else:
    #         raise ValueError(f"Unknown model '{spec.model}'")

    def _fit_and_blockdrop(
        Xmat: np.ndarray,
        yvec: np.ndarray,
        mask: Optional[np.ndarray] = None,
    ) -> Tuple[float, np.ndarray, float, Dict[str, float]]:
        if mask is not None:
            mask = np.asarray(mask, dtype=bool)
            if mask.ndim != 1 or mask.shape[0] != X.shape[0]:
                raise ValueError(f"mask shape {mask.shape} != (T,) with T={X.shape[0]}")
            Xmat = Xmat[mask]
            yvec = yvec[mask]
        if spec.model == "poisson":
            return fit_unit_poisson_cv(Xmat, yvec, alpha=spec.alpha, n_splits=spec.n_splits, seed=spec.random_seed)
        elif spec.model == "ridge":
            # common choice: z-score y for ridge
            return fit_unit_ridge_cv(Xmat, _safe_zscore(yvec), alpha=spec.alpha, n_splits=spec.n_splits, seed=spec.random_seed)
        else:
            raise ValueError(f"Unknown model '{spec.model}'")
        
    for n in range(N):
        y = Yg[:, n].astype(np.float32)
        # if constant / all-zero, skip
        if np.all(~np.isfinite(y)) or np.nanstd(y) < 1e-8:
            r2_full[n] = np.nan
            coef_full[n, :] = np.nan
            intercept_full[n] = np.nan
            for bn in delta_r2.keys():
                delta_r2[bn][n] = np.nan
            continue

        r2, coef, icpt = _fit_and_blockdrop(Xg, y, mask=fit_mask)
        r2_full[n] = r2
        coef_full[n, :] = coef
        intercept_full[n] = icpt

        # block-drop: refit without the block and compute delta
        for bn, Xred in reduced_X.items():
            r2_red, _, _ = _fit_and_blockdrop(Xred, y, mask=fit_mask)
            delta_r2[bn][n] = r2_full[n] - r2_red

    # ---- Save
    meta_out = {
        "bin_size_s": float(bin_size_s),
        "T_bins": int(T),
        "T_used": int(Tg),
        "N_units": int(N),
        "model": spec.model,
        "alpha": float(spec.alpha),
        "cv_n_splits": int(spec.n_splits),
        "drift_condition": drift_condition,
        "drift_pcs_k": int(spec.drift_pcs_k),
        "stim_f_hz": float(spec.stim_f_hz),
        "stim_impulse_kernel_s": float(spec.stim_impulse_kernel_s),
        "include_phi2": bool(spec.include_phi2),
        "behavior_keys": list(spec.behavior_keys),
        "design_names": names,
        "block_slices": {k: list(v) for k, v in slices.items()},
    }

    with h5py.File(out_h5_path, "w") as f:
        f.attrs["meta_json"] = json.dumps(meta_out)

        # design
        gD = f.create_group("design")
        gD.create_dataset("X_names", data=np.array(names, dtype="S"))
        for bn in ["A_drift", "B_behavior", "C_stim", "D_phase"]:
            g = gD.create_group(bn)
            g.create_dataset("names", data=np.array(names_blocks[bn], dtype="S"))
            g.create_dataset("slice", data=np.array(slices[bn], dtype=np.int64))

        # masks / time
        f.create_dataset("time/bin_centers_s", data=bin_centers_s.astype(np.float32))
        f.create_dataset("time/good_row_mask", data=good_row.astype(np.uint8))
        f.create_dataset("states/mask_target", data=mask_tgt.astype(np.uint8))

        # fits
        gF = f.create_group("fit")
        gF.create_dataset("coef_full", data=coef_full)          # (N, P)
        gF.create_dataset("intercept_full", data=intercept_full)  # (N,)
        gF.create_dataset("r2_full_cv", data=r2_full)           # (N,)

        gDrop = f.create_group("fit/block_drop")
        for bn, arr in delta_r2.items():
            gDrop.create_dataset(f"delta_r2_{bn}", data=arr)

        # optional unit meta placeholder
        if unit_meta is not None:
            f.create_dataset("unit_meta", data=np.asarray(unit_meta))

    #print(f"[cell_glm] wrote {out_h5_path}")
