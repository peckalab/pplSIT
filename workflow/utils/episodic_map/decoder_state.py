import json
import numpy as np
import h5py

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix, accuracy_score


# ----------------------------
# Utilities
# ----------------------------
def as_sorted_unique_int(arr):
    arr = np.asarray(arr, dtype=int)
    arr = arr[np.isfinite(arr)]
    arr = np.unique(arr)
    arr = arr[arr >= 0]
    return np.sort(arr)


def make_block_folds(n_samples, n_splits=5, min_test_frac=0.15):
    """
    Blocked CV: split contiguous samples into n_splits blocks.
    Each fold uses one block as test, rest as train.
    """
    if n_splits < 2:
        raise ValueError("n_splits must be >= 2")
    if n_samples < 10:
        raise ValueError("Too few samples for CV")

    # split indices into contiguous blocks
    edges = np.linspace(0, n_samples, n_splits + 1).round().astype(int)
    folds = []
    for k in range(n_splits):
        te = np.arange(edges[k], edges[k + 1])
        tr = np.setdiff1d(np.arange(n_samples), te, assume_unique=False)
        if len(te) / n_samples < min_test_frac and n_splits > 2:
            # if last block too small due to rounding, merge with previous
            continue
        folds.append((tr, te))
    if len(folds) < 2:
        # fallback: even/odd split
        idx = np.arange(n_samples)
        folds = [(idx[idx % 2 == 0], idx[idx % 2 == 1])]
    return folds


def zscore_train_apply(X_train, X_test, eps=1e-6):
    """
    Z-score per feature using train stats only.
    X_* are (n_samples, n_features).
    """
    mu = X_train.mean(axis=0, keepdims=True)
    sd = X_train.std(axis=0, keepdims=True)
    sd = np.maximum(sd, eps)
    return (X_train - mu) / sd, (X_test - mu) / sd, mu.squeeze(), sd.squeeze()


def balanced_subsample_indices(y, rng, n_per_class=None):
    """
    Returns indices to make class counts equal.
    If n_per_class is None: uses min class count.
    """
    y = np.asarray(y, dtype=int)
    classes = np.unique(y)
    idxs_by_class = {c: np.where(y == c)[0] for c in classes}
    counts = {c: len(idxs_by_class[c]) for c in classes}
    if n_per_class is None:
        n_per_class = min(counts.values())
    keep = []
    for c in classes:
        if counts[c] < n_per_class:
            raise ValueError(f"Class {c} has only {counts[c]} samples < n_per_class={n_per_class}")
        keep.append(rng.choice(idxs_by_class[c], size=n_per_class, replace=False))
    keep = np.concatenate(keep)
    keep = rng.permutation(keep)
    return keep, n_per_class


def fit_predict_fold(X, y, tr, te, C=1.0, max_iter=2000):
    """
    Multinomial logistic regression.
    """
    Xtr, Xte = X[tr], X[te]
    ytr, yte = y[tr], y[te]

    Xtr_z, Xte_z, mu, sd = zscore_train_apply(Xtr, Xte)

    clf = LogisticRegression(
        #multi_class="multinomial",
        solver="lbfgs",
        C=C,
        max_iter=max_iter,
        #n_jobs=1,
    )
    clf.fit(Xtr_z, ytr)
    yhat = clf.predict(Xte_z)
    return yte, yhat, mu, sd

