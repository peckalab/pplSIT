import numpy as np
import h5py
from sklearn.decomposition import PCA
import json


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


