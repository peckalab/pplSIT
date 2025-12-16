import os, sys
import numpy as np

from sklearn.neighbors import NearestNeighbors


STATE_KEYS = [
    "idxs_tgt_sta_succ",
    "idxs_bgr_sta",
    "idxs_bgr_run",
    "idxs_sil_sta",
    "idxs_sil_run",
]

STATE_NAMES  = ["tgt", "bgr_sta", "bgr_run", "sil_sta", "sil_run"]
STATE_LABELS = ['Target', 'BGR/sta', 'BGR/run', 'No stim./sta', 'No stim./run']
STATE_COLORS = ['tab:orange', 'mediumblue', 'deepskyblue', 'black', 'grey']


def build_state_labels(t_bins: int, state_idxs: dict, state_keys=STATE_KEYS):
    """
    labels: int array shape (t_bins,), -1 = unlabeled, else 0..4 in STATE_KEYS order.
    """
    labels = np.full(t_bins, -1, dtype=np.int16)

    # Assign labels (later keys overwrite earlier ones if overlap exists)
    for si, k in enumerate(state_keys):
        if k not in state_idxs:
            raise KeyError(f"Missing key in state_idxs: {k}")
        idx = np.asarray(state_idxs[k], dtype=np.int64)
        if idx.size == 0:
            continue
        if idx.min() < 0 or idx.max() >= t_bins:
            raise ValueError(f"Indices for {k} out of range [0, {t_bins - 1}]")
        labels[idx] = si

    # Optional overlap warning (pairwise)
    overlaps = []
    for i in range(len(state_keys)):
        a = set(map(int, np.asarray(state_idxs[state_keys[i]], dtype=np.int64)))
        for j in range(i + 1, len(state_keys)):
            b = set(map(int, np.asarray(state_idxs[state_keys[j]], dtype=np.int64)))
            inter = a.intersection(b)
            if inter:
                overlaps.append((state_keys[i], state_keys[j], len(inter)))

    return labels, overlaps


# -----------------------------
# Preprocessing
# -----------------------------
def preprocess_counts(
    X_counts: np.ndarray,
    X_is_neurons_by_time: bool = True,
    transform: str = "sqrt",
    zscore: bool = True,
):
    """
    Input:
      X_counts: spike counts, shape (n_neurons, t_bins) if X_is_neurons_by_time else (t_bins, n_neurons)

    Output:
      Xz: (t_bins, n_neurons) float32
      stats: dict with per-neuron mean/std (on transformed data)
    """
    if X_is_neurons_by_time:
        X = X_counts.T
    else:
        X = X_counts.copy()

    X = X.astype(np.float64)

    if transform == "sqrt":
        X = np.sqrt(X)
    elif transform == "log1p":
        X = np.log1p(X)
    elif transform is None:
        pass
    else:
        raise ValueError("transform must be 'sqrt', 'log1p', or None")

    stats = {}
    if zscore:
        mu = X.mean(axis=0, keepdims=True)
        sd = X.std(axis=0, keepdims=True)
        sd[sd == 0] = 1.0
        X = (X - mu) / sd
        stats["mu"] = mu.squeeze().astype(np.float32)
        stats["sd"] = sd.squeeze().astype(np.float32)

    # Drop non-finite rows (rare but protects embeddings)
    finite_mask = np.isfinite(X).all(axis=1)
    return X[finite_mask].astype(np.float32), stats, finite_mask


# -----------------------------
# kNN purity + nulls
# -----------------------------
def knn_purity(Y, labels, target_label: int, k: int, exclude_self: bool = True) -> float:
    Y = np.asarray(Y, dtype=np.float64)
    labels = np.asarray(labels)

    idx_t = np.where(labels == target_label)[0]
    if idx_t.size == 0:
        return np.nan

    n_neighbors = k + 1 if exclude_self else k
    nn = NearestNeighbors(n_neighbors=n_neighbors, metric="euclidean")
    nn.fit(Y)
    neigh = nn.kneighbors(Y[idx_t], return_distance=False)

    if exclude_self:
        # Remove self neighbor robustly
        rows = []
        for qi, row in zip(idx_t, neigh):
            row = row[row != qi]
            rows.append(row[:k])
        neigh = np.vstack(rows)

    return float((labels[neigh] == target_label).mean())


def null_circular_shift(
    Y,
    labels,
    target_label: int,
    k: int,
    n_perm: int,
    seed: int = 0,
    min_shift: int = 10,
):
    """
    Circularly shift labels to preserve temporal contiguity structure.
    Returns array of null purity values (len n_perm).
    """
    rng = np.random.default_rng(seed)
    labels = np.asarray(labels)
    n = len(labels)

    null_vals = np.empty(n_perm, dtype=np.float32)
    for i in range(n_perm):
        shift = int(rng.integers(min_shift, n - min_shift))
        lab_shift = np.roll(labels, shift)
        null_vals[i] = knn_purity(Y, lab_shift, target_label=target_label, k=k)
    return null_vals


def summarize_vs_null(obs: float, null_vals: np.ndarray):
    null_vals = np.asarray(null_vals, dtype=np.float64)
    mu = float(np.nanmean(null_vals))
    sd = float(np.nanstd(null_vals, ddof=1)) if np.isfinite(null_vals).sum() > 1 else np.nan
    # One-sided p-value with add-one smoothing
    p = float((np.sum(null_vals >= obs) + 1) / (len(null_vals) + 1))
    z = float((obs - mu) / sd) if (sd is not None and sd > 0) else np.nan
    return obs, mu, sd, z, p


def compute_knn_stats(
    Y2,
    labels,
    target_label: int,
    k_list=(15, 30, 50),
    n_perm: int = 1000,
    seed: int = 0,
):
    """
    Returns dict keyed by k with (obs, null_mean, null_std, z, p, null_vals)
    """
    out = {}
    for k in k_list:
        obs = knn_purity(Y2, labels, target_label=target_label, k=k)
        null_vals = null_circular_shift(
            Y2, labels, target_label=target_label, k=k, n_perm=n_perm, seed=seed
        )
        obs, mu, sd, z, p = summarize_vs_null(obs, null_vals)
        out[int(k)] = {
            "obs": float(obs),
            "null_mean": float(mu),
            "null_std": float(sd),
            "z": float(z),
            "p_one_sided": float(p),
            "null_vals": null_vals.astype(np.float32),
        }
    return out
