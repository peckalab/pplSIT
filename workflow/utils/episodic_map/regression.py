import json
import numpy as np
import h5py
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

from sklearn.linear_model import LinearRegression, Ridge
from sklearn.metrics import r2_score


# ----------------------------
# Utilities
# ----------------------------

def upper_tri_indices(n: int) -> Tuple[np.ndarray, np.ndarray]:
    """Return i,j indices for upper triangle (i<j)."""
    return np.triu_indices(n, k=1)

def upper_tri_vec(M: np.ndarray) -> np.ndarray:
    """Vectorize upper triangle (i<j)."""
    i, j = upper_tri_indices(M.shape[0])
    return M[i, j]

def zscore_cols(X: np.ndarray, eps: float = 1e-8) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Z-score each column. Returns Xz, mu, sd."""
    mu = np.nanmean(X, axis=0)
    sd = np.nanstd(X, axis=0)
    sd = np.where(sd < eps, 1.0, sd)
    Xz = (X - mu) / sd
    return Xz, mu, sd

def safe_mask(*vecs: np.ndarray) -> np.ndarray:
    """Mask for finite entries across all provided vectors."""
    m = np.ones_like(vecs[0], dtype=bool)
    for v in vecs:
        m &= np.isfinite(v)
    return m

def permute_matrix_by_episode(M: np.ndarray, perm: np.ndarray) -> np.ndarray:
    """Apply the same episode permutation to rows and columns."""
    return M[np.ix_(perm, perm)]

def partial_r2_from_full_and_reduced(y: np.ndarray, yhat_full: np.ndarray, yhat_red: np.ndarray) -> float:
    """
    Partial R^2 for a predictor group:
      (SSE_reduced - SSE_full) / SSE_reduced
    """
    sse_full = np.sum((y - yhat_full) ** 2)
    sse_red  = np.sum((y - yhat_red) ** 2)
    if sse_red <= 1e-12:
        return np.nan
    return float((sse_red - sse_full) / sse_red)

def _slice_square(M: np.ndarray, n: int) -> np.ndarray:
    M = np.asarray(M)
    if M.ndim != 2 or M.shape[0] != M.shape[1]:
        raise ValueError(f"Expected square matrix, got {M.shape}")
    return M[:n, :n]


@dataclass
class RegressionSpec:
    # Which similarity to use as y
    y_key: str = "S_latent"     # or "S_pv"
    # Predictors to include
    # Added:
    #   - dcenter_m: |dist_to_center_i - dist_to_center_j|
    #   - dhd_rel_center_rad: circular distance between HD relative to center
    x_keys: Tuple[str, ...] = (
        "dt_s",
        "dspace_m",
        "dcenter_m",
        "dصدhd_rad" if False else "dhd_rad",  # no-op; remove if you want (keeps diff minimal)
        "dhd_rel_center_rad",
        "dturn",
        "dstill",
    )
    # Model type
    model: str = "ridge"        # "ols" or "ridge"
    ridge_alpha: float = 1.0
    # Null settings
    n_perm_episode_shuffle: int = 200
    n_perm_y_shuffle: int = 200
    n_perm_circshift: int = 0   # optional
    random_seed: int = 0


# ----------------------------
# HDF5 I/O
# ----------------------------

def load_episode_similarity_session(h5_path: str) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
    """
    Loads:
      - similarity matrices from episode_similarity/similarity/
      - covariate matrices from episode_similarity/covariates/
    Returns:
      sim: dict of similarity matrices
      cov: dict of covariate matrices
    """
    sim = {}
    cov = {}
    with h5py.File(h5_path, "r") as f:
        g_sim = f["episode_similarity/similarity"]
        for k in g_sim.keys():
            sim[k] = np.asarray(g_sim[k][...])

        g_cov = f["episode_similarity/covariates"]
        for k in g_cov.keys():
            cov[k] = np.asarray(g_cov[k][...])

    return sim, cov


def write_regression_results_h5(
    out_h5: str,
    *,
    meta: Dict,
    betas: np.ndarray,
    beta_names: List[str],
    intercept: float,
    r2_full: float,
    partial_r2: Dict[str, float],
    n_pairs: int,
    y_key: str,
    x_keys: List[str],
    nulls: Dict[str, Dict[str, np.ndarray]],
):
    """
    Saves regression summary + null distributions.

    nulls format:
      nulls["episode_shuffle"]["r2"] -> (n_perm,)
      nulls["episode_shuffle"]["betas"] -> (n_perm, n_x)
      etc.
    """
    with h5py.File(out_h5, "w") as f:
        f.attrs["meta_json"] = json.dumps(meta)

        g = f.create_group("regression")
        g.attrs["y_key"] = y_key
        g.attrs["x_keys_json"] = json.dumps(list(x_keys))
        g.create_dataset("beta_names", data=np.array(beta_names, dtype="S"), compression="gzip")

        g.create_dataset("betas", data=betas.astype(np.float32), compression="gzip")
        g.create_dataset("intercept", data=np.float32(intercept))  # scalar: no compression
        g.create_dataset("r2_full", data=np.float32(r2_full))
        g.create_dataset("n_pairs", data=np.int64(n_pairs))

        g_pr2 = g.create_group("partial_r2")
        for k, v in partial_r2.items():
            g_pr2.create_dataset(k, data=np.float32(v))

        g_null = f.create_group("nulls")
        for null_name, dd in nulls.items():
            gn = g_null.create_group(null_name)
            for k, arr in dd.items():
                if np.ndim(arr) == 0:
                    gn.create_dataset(k, data=np.asarray(arr))
                else:
                    gn.create_dataset(k, data=np.asarray(arr), compression="gzip")


# ----------------------------
# Core regression
# ----------------------------

def fit_model(X: np.ndarray, y: np.ndarray, spec: RegressionSpec):
    if spec.model == "ols":
        model = LinearRegression()
    elif spec.model == "ridge":
        model = Ridge(alpha=spec.ridge_alpha, fit_intercept=True)
    else:
        raise ValueError(f"Unknown model: {spec.model}")
    model.fit(X, y)
    yhat = model.predict(X)
    return model, yhat

def run_session_regression(
    h5_path: str,
    out_h5: str,
    spec: RegressionSpec,
):
    rng = np.random.default_rng(spec.random_seed)

    sim, cov = load_episode_similarity_session(h5_path)
    if spec.y_key not in sim:
        raise KeyError(f"Missing y_key {spec.y_key} in episode_similarity/similarity. Found: {list(sim.keys())}")

    # Load y matrix and predictors
    Y = np.asarray(sim[spec.y_key])
    if Y.ndim != 2 or Y.shape[0] != Y.shape[1]:
        raise ValueError(f"Similarity matrix must be square. Got {Y.shape}")

    X_mats = {}
    for k in spec.x_keys:
        if k not in cov:
            raise KeyError(f"Missing covariate {k} in episode_similarity/covariates. Found: {list(cov.keys())}")
        M = np.asarray(cov[k])
        if M.ndim != 2 or M.shape[0] != M.shape[1]:
            raise ValueError(f"Covariate {k} must be square. Got {M.shape}")
        X_mats[k] = M

    # --- NEW: reconcile episode counts across Y and covariates ---
    sizes = [Y.shape[0]] + [X_mats[k].shape[0] for k in spec.x_keys]
    n_ep = int(np.min(sizes))

    if any(s != n_ep for s in sizes):
        print(
            f"[run_session_regression] WARNING size mismatch in {h5_path}. "
            f"Using n_ep={n_ep}. Y:{Y.shape} " +
            " ".join([f"{k}:{X_mats[k].shape}" for k in spec.x_keys])
        )

    Y = _slice_square(Y, n_ep)
    for k in list(X_mats.keys()):
        X_mats[k] = _slice_square(X_mats[k], n_ep)

    # Vectorize pairs
    y_vec = upper_tri_vec(Y)
    x_vecs = {k: upper_tri_vec(M) for k, M in X_mats.items()}

    # Valid mask across all vectors
    m = safe_mask(y_vec, *[x_vecs[k] for k in spec.x_keys])
    y = y_vec[m].astype(float)

    X = np.column_stack([x_vecs[k][m].astype(float) for k in spec.x_keys])

    # Z-score predictors (recommended)
    Xz, x_mu, x_sd = zscore_cols(X)

    # Fit full model
    model_full, yhat_full = fit_model(Xz, y, spec)
    r2_full = float(r2_score(y, yhat_full))

    # Standardized betas (since X is z-scored)
    betas = np.asarray(model_full.coef_, dtype=float)
    intercept = float(model_full.intercept_)

    # Partial R^2 per predictor (leave-one-out reduced model)
    partial_r2 = {}
    for j, key in enumerate(spec.x_keys):
        cols = [c for c in range(Xz.shape[1]) if c != j]
        if len(cols) == 0:
            partial_r2[key] = np.nan
            continue
        _, yhat_red = fit_model(Xz[:, cols], y, spec)
        partial_r2[key] = partial_r2_from_full_and_reduced(y, yhat_full, yhat_red)

    # -------- Nulls --------
    nulls = {}

    # Precompute covariate vectors once (fixed X)
    x_vecs = {k: upper_tri_vec(M) for k, M in X_mats.items()}

    # (1) Episode-label shuffle NULL (FIXED):
    # permute episode indices in Y ONLY; keep X fixed.
    if spec.n_perm_episode_shuffle > 0:
        r2s = np.zeros(spec.n_perm_episode_shuffle, dtype=np.float32)
        bss = np.zeros((spec.n_perm_episode_shuffle, len(spec.x_keys)), dtype=np.float32)

        for p in range(spec.n_perm_episode_shuffle):
            perm = rng.permutation(n_ep)

            # Permute only Y
            Yp = permute_matrix_by_episode(Y, perm)
            yv = upper_tri_vec(Yp)

            # Use original covariates (unpermuted)
            mp = safe_mask(yv, *[x_vecs[k] for k in spec.x_keys])

            yy = yv[mp].astype(float)
            XX = np.column_stack([x_vecs[k][mp].astype(float) for k in spec.x_keys])

            XXz, _, _ = zscore_cols(XX)

            mod, yhat = fit_model(XXz, yy, spec)
            r2s[p] = r2_score(yy, yhat)
            bss[p] = mod.coef_.astype(np.float32)

        nulls["episode_shuffle"] = {"r2": r2s, "betas": bss}


    # (2) y-shuffle across pairs: already correct
    if spec.n_perm_y_shuffle > 0:
        r2s = np.zeros(spec.n_perm_y_shuffle, dtype=np.float32)
        for p in range(spec.n_perm_y_shuffle):
            yy = rng.permutation(y)
            mod, yhat = fit_model(Xz, yy, spec)
            r2s[p] = r2_score(yy, yhat)
        nulls["y_shuffle"] = {"r2": r2s}


    # (3) Circular shift NULL (FIXED):
    # circularly shift episode indices in Y ONLY; keep X fixed.
    if spec.n_perm_circshift > 0:
        r2s = np.zeros(spec.n_perm_circshift, dtype=np.float32)
        bss = np.zeros((spec.n_perm_circshift, len(spec.x_keys)), dtype=np.float32)
        shifts = rng.integers(low=1, high=n_ep, size=spec.n_perm_circshift)

        for p, sh in enumerate(shifts):
            perm = (np.arange(n_ep) + int(sh)) % n_ep

            # Permute only Y
            Yp = permute_matrix_by_episode(Y, perm)
            yv = upper_tri_vec(Yp)

            # Use original covariates (unpermuted)
            mp = safe_mask(yv, *[x_vecs[k] for k in spec.x_keys])

            yy = yv[mp].astype(float)
            XX = np.column_stack([x_vecs[k][mp].astype(float) for k in spec.x_keys])

            XXz, _, _ = zscore_cols(XX)

            mod, yhat = fit_model(XXz, yy, spec)
            r2s[p] = r2_score(yy, yhat)
            bss[p] = mod.coef_.astype(np.float32)

        nulls["circular_shift"] = {"r2": r2s, "betas": bss, "shift": shifts.astype(np.int64)}

    # Meta
    meta = dict(
        source_h5=h5_path,
        y_key=spec.y_key,
        x_keys=list(spec.x_keys),
        model=spec.model,
        ridge_alpha=float(spec.ridge_alpha),
        n_ep=int(n_ep),
        n_pairs=int(len(y)),
        x_mu=x_mu.tolist(),
        x_sd=x_sd.tolist(),
        random_seed=int(spec.random_seed),
        n_perm_episode_shuffle=int(spec.n_perm_episode_shuffle),
        n_perm_y_shuffle=int(spec.n_perm_y_shuffle),
        n_perm_circshift=int(spec.n_perm_circshift),
    )

    write_regression_results_h5(
        out_h5,
        meta=meta,
        betas=betas,
        beta_names=list(spec.x_keys),
        intercept=intercept,
        r2_full=r2_full,
        partial_r2=partial_r2,
        n_pairs=len(y),
        y_key=spec.y_key,
        x_keys=list(spec.x_keys),
        nulls=nulls
    )

    return out_h5
