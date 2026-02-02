import numpy as np
import h5py


def _gaussian_smooth_1d(x: np.ndarray, sigma_bins: float) -> np.ndarray:
    """Gaussian smoothing with reflect padding. sigma in bins."""
    if sigma_bins is None or sigma_bins <= 0:
        return x.astype(float, copy=True)

    # kernel length ~ +/- 4 sigma
    half = int(np.ceil(4 * sigma_bins))
    kx = np.arange(-half, half + 1)
    k = np.exp(-0.5 * (kx / sigma_bins) ** 2)
    k /= k.sum()

    # reflect-pad and convolve
    xp = np.pad(x, (half, half), mode="reflect")
    y = np.convolve(xp, k, mode="valid")
    return y

def _robust_zscore(x: np.ndarray, eps: float = 1e-8) -> np.ndarray:
    """Robust z-score using median and MAD (scaled)."""
    med = np.nanmedian(x)
    mad = np.nanmedian(np.abs(x - med))
    robust_std = 1.4826 * mad  # MAD -> std for Gaussian
    return (x - med) / (robust_std + eps)

def _enforce_min_duration(mask: np.ndarray, min_len: int) -> np.ndarray:
    """Keep only True runs of length >= min_len."""
    if min_len <= 1:
        return mask

    m = mask.astype(bool)
    out = np.zeros_like(m, dtype=bool)

    # run-length encoding
    idx = np.flatnonzero(np.diff(np.r_[False, m, False]))
    # idx comes in pairs: start, end
    starts = idx[0::2]
    ends = idx[1::2]
    for s, e in zip(starts, ends):
        if (e - s) >= min_len:
            out[s:e] = True
    return out

def build_gate_engagement_mask(
    gate_counts_nt: np.ndarray,
    *,
    dt_s: float = 0.05,
    smooth_sigma_bins: float = 2.0,     # ~100 ms at 50 ms bins
    thr_quantile: float = 0.95,         # robust threshold on ensemble score
    hysteresis_delta: float = 0.5,      # in z units of ensemble score
    min_on_bins: int = 2,               # minimum ON duration (e.g. 2 bins = 100 ms)
    eps: float = 1e-8,
) -> dict:
    """
    Variant A: "Ensemble engagement" mask for gate cells, assuming all ramp UP.

    Input
    -----
    gate_counts_nt : array (n_units, n_timebins)
        Binned spike counts (or rates) for gate units. Can be (1, T).

    Steps
    -----
    1) Smooth each unit (Gaussian).
    2) Robust-normalize each unit (median/MAD z-score across full session).
    3) Ensemble score = mean across units.
    4) Threshold using quantile of the ensemble score distribution.
    5) Hysteresis + minimum ON duration.

    Returns
    -------
    dict with:
      - "mask": bool array (T,)
      - "score": float array (T,) ensemble engagement score (z-ish units)
      - "thr_on", "thr_off": thresholds used
    """
    X = np.asarray(gate_counts_nt, dtype=float)
    if X.ndim != 2:
        raise ValueError(f"gate_counts_nt must be 2D (n_units, T). Got shape {X.shape}")
    n_units, T = X.shape
    if T < 5:
        raise ValueError("Too few time bins to build a stable mask.")

    # 1) Smooth each unit
    Xs = np.zeros_like(X, dtype=float)
    for i in range(n_units):
        Xs[i] = _gaussian_smooth_1d(X[i], sigma_bins=smooth_sigma_bins)

    # 2) Robust z-score per unit across the whole session
    Xz = np.zeros_like(Xs, dtype=float)
    for i in range(n_units):
        Xz[i] = _robust_zscore(Xs[i], eps=eps)

    # 3) Ensemble engagement score (assume coherent ramp UP)
    score = np.nanmean(Xz, axis=0)

    # 4) Robust thresholds from score distribution
    thr_on = float(np.nanquantile(score, thr_quantile))
    thr_off = thr_on - float(hysteresis_delta)

    # 5) Apply hysteresis to get a stable mask
    mask = np.zeros(T, dtype=bool)
    state = False
    for t in range(T):
        if not state:
            if score[t] >= thr_on:
                state = True
        else:
            if score[t] <= thr_off:
                state = False
        mask[t] = state

    # Enforce minimum ON duration
    mask = _enforce_min_duration(mask, min_on_bins)

    return {"mask": mask, "score": score, "thr_on": thr_on, "thr_off": thr_off}


def save_gate_engagement_to_h5(
    h5_path: str,
    h5_group: str,
    out: dict,
    *,
    overwrite: bool = True,
):
    """
    Save gate engagement mask + score to HDF5.

    Parameters
    ----------
    h5_path : str
        Path to HDF5 file (created if not exists)
    h5_group : str
        Group path inside HDF5, e.g.:
            "/gate_masks/variantA/session"
    out : dict
        Output of build_gate_engagement_mask()
    overwrite : bool
        If True, overwrite existing datasets in the group
    """
    with h5py.File(h5_path, "a") as f:
        # create or open group
        if h5_group in f:
            g = f[h5_group]
            if overwrite:
                for k in list(g.keys()):
                    del g[k]
        else:
            g = f.create_group(h5_group)

        # save arrays
        g.create_dataset(
            "mask",
            data=out["mask"].astype(np.uint8),  # compact, boolean semantics
            compression="gzip",
            compression_opts=4,
        )
        g.create_dataset(
            "score",
            data=out["score"].astype(np.float32),
            compression="gzip",
            compression_opts=4,
        )

        # save thresholds as attributes
        g.attrs["thr_on"] = float(out["thr_on"])
        g.attrs["thr_off"] = float(out["thr_off"])
