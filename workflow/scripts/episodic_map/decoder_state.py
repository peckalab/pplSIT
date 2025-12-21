import os, sys
import h5py
import json
import numpy as np


# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.episodic_map.decoder_state import *

# ----------------------------
# Main single-session routine
# ----------------------------
def run_stationary_state_decoder_single_session(
    X_counts,
    state_idxs,
    t_edges=None,
    out_h5_path="state_decoder.h5",
    keys=("idxs_sil_sta", "idxs_bgr_sta", "idxs_tgt_sta_succ"),
    label_names=("sil_sta", "bgr_sta", "tgt_sta"),
    n_splits=5,
    C=1.0,
    max_iter=2000,
    do_balance=True,
    balance_n_per_class=None,   # None -> min count
    do_label_shuffle=True,
    n_label_shuffles=50,
    rng_seed=0,
):
    """
    X_counts: (n_units, n_bins) spike counts at 250ms.
    state_idxs: dict key->array of bin indices at same sampling.
    t_edges: optional, (n_bins+1,)
    """

    rng = np.random.default_rng(rng_seed)

    X_counts = np.asarray(X_counts)
    if X_counts.ndim != 2:
        raise ValueError("X_counts must be (n_units, n_bins)")
    n_units, n_bins = X_counts.shape

    # Build samples: X (n_samples, n_units), y (n_samples,)
    all_bins = []
    all_labels = []

    # map labels 0..K-1
    for lab, k in enumerate(keys):
        if k not in state_idxs:
            raise KeyError(f"Missing state key in state_idxs: {k}")
        idxs = as_sorted_unique_int(state_idxs[k])
        idxs = idxs[idxs < n_bins]
        if len(idxs) == 0:
            raise ValueError(f"State {k} has 0 bins after filtering")
        all_bins.append(idxs)
        all_labels.append(np.full(len(idxs), lab, dtype=int))

    bins = np.concatenate(all_bins)
    y = np.concatenate(all_labels)

    # Sort by time (important for blocked CV)
    order = np.argsort(bins)
    bins = bins[order]
    y = y[order]

    X = X_counts[:, bins].T.astype(np.float32)  # (n_samples, n_units)

    # Balance (global, before CV) to keep folds comparable
    balance_info = {}
    if do_balance:
        keep, n_per_class = balanced_subsample_indices(y, rng, n_per_class=balance_n_per_class)
        X = X[keep]
        y = y[keep]
        bins = bins[keep]
        # re-sort by time after balancing
        o2 = np.argsort(bins)
        X, y, bins = X[o2], y[o2], bins[o2]
        balance_info = {"enabled": True, "n_per_class": int(n_per_class)}
    else:
        balance_info = {"enabled": False}

    n_samples = X.shape[0]
    folds = make_block_folds(n_samples, n_splits=n_splits)

    n_classes = len(keys)
    cm_sum = np.zeros((n_classes, n_classes), dtype=np.float64)
    accs = []

    # Store mu/sd per fold? (optional; can be large). We'll store only fold accuracies + cm.
    for tr, te in folds:
        yte, yhat, _, _ = fit_predict_fold(X, y, tr, te, C=C, max_iter=max_iter)
        cm_sum += confusion_matrix(yte, yhat, labels=np.arange(n_classes))
        accs.append(accuracy_score(yte, yhat))

    cm_mean = cm_sum / cm_sum.sum(axis=1, keepdims=True).clip(min=1.0)  # row-normalized
    acc_mean = float(np.mean(accs))
    acc_std = float(np.std(accs, ddof=1)) if len(accs) > 1 else float("nan")

    # Label-shuffle null
    shuf_accs = []
    shuf_cm_sum = np.zeros((n_classes, n_classes), dtype=np.float64)

    if do_label_shuffle:
        # keep the same folds; shuffle y across samples (destroys state structure)
        for s in range(n_label_shuffles):
            y_shuf = rng.permutation(y)

            cm_s = np.zeros((n_classes, n_classes), dtype=np.float64)
            acc_s = []

            for tr, te in folds:
                yte, yhat, _, _ = fit_predict_fold(X, y_shuf, tr, te, C=C, max_iter=max_iter)
                cm_s += confusion_matrix(yte, yhat, labels=np.arange(n_classes))
                acc_s.append(accuracy_score(yte, yhat))

            shuf_cm_sum += cm_s / cm_s.sum(axis=1, keepdims=True).clip(min=1.0)
            shuf_accs.append(np.mean(acc_s))

        shuf_acc_mean = float(np.mean(shuf_accs))
        shuf_acc_std = float(np.std(shuf_accs, ddof=1)) if len(shuf_accs) > 1 else float("nan")
        shuf_cm_mean = (shuf_cm_sum / len(shuf_accs)).astype(np.float32)
    else:
        shuf_acc_mean = float("nan")
        shuf_acc_std = float("nan")
        shuf_cm_mean = np.full((n_classes, n_classes), np.nan, dtype=np.float32)

    # Store to HDF5
    meta = {
        "label_names": list(label_names),
        "keys": list(keys),
        "n_units": int(n_units),
        "n_bins_total": int(n_bins),
        "n_samples_used": int(n_samples),
        "cv": {
            "scheme": "blocked",
            "n_splits": int(n_splits),
            "n_folds_effective": int(len(folds)),
        },
        "classifier": {
            "type": "multinomial_logreg",
            "C": float(C),
            "max_iter": int(max_iter),
        },
        "balance": balance_info,
        "label_shuffle": {
            "enabled": bool(do_label_shuffle),
            "n_shuffles": int(n_label_shuffles) if do_label_shuffle else 0,
        },
        "rng_seed": int(rng_seed),
    }

    with h5py.File(out_h5_path, "w") as h5:
        h5.attrs["meta_json"] = json.dumps(meta)

        h5.create_dataset("data/bins_used", data=bins.astype(np.int64))
        h5.create_dataset("data/y_labels", data=y.astype(np.int16))
        h5.create_dataset("results/acc_folds", data=np.array(accs, dtype=np.float32))
        h5.create_dataset("results/acc_mean", data=np.array(acc_mean, dtype=np.float32))
        h5.create_dataset("results/acc_std", data=np.array(acc_std, dtype=np.float32))
        h5.create_dataset("results/confusion_row_norm", data=cm_mean.astype(np.float32))

        h5.create_dataset("null/acc_shuffles", data=np.array(shuf_accs, dtype=np.float32))
        h5.create_dataset("null/acc_shuf_mean", data=np.array(shuf_acc_mean, dtype=np.float32))
        h5.create_dataset("null/acc_shuf_std", data=np.array(shuf_acc_std, dtype=np.float32))
        h5.create_dataset("null/confusion_row_norm_mean", data=shuf_cm_mean)

        if t_edges is not None:
            t_edges = np.asarray(t_edges)
            if t_edges.ndim == 1 and len(t_edges) == n_bins + 1:
                h5.create_dataset("data/t_edges", data=t_edges.astype(np.float64))

    return {
        "out_h5_path": out_h5_path,
        "acc_mean": acc_mean,
        "acc_std": acc_std,
        "acc_shuf_mean": shuf_acc_mean,
        "acc_shuf_std": shuf_acc_std,
    }

#cfg = snakemake.config['trajectories']

STATE_KEYS = [
    "idxs_tgt_sta_succ",
    "idxs_bgr_sta",
    "idxs_sil_sta",
]

# -----------------------------
# Reading datasets
# -----------------------------
meta_file = snakemake.input[0]
segm_file = snakemake.input[1]
actm_file = snakemake.input[2]

with h5py.File(meta_file, 'r') as f:
    tl = np.array(f['processed']['timeline'])
    tgt_mx = np.array(f['processed']['target_matrix'])
    sound_events = np.array(f['processed']['sound_events'])

with h5py.File(actm_file, 'r') as f:
    X_counts = np.array(f['mx_250ms']['mx'])
    t_edges  = np.array(f['mx_250ms']['bins'])
    
#X_counts_50ms = gaussian_filter1d(X_counts_50ms, sigma=5, axis=0, mode="nearest")

with h5py.File(segm_file, 'r') as f:
    state_idxs = {}
    for idxs_name in STATE_KEYS:
        state_idxs[idxs_name] = np.array(f[idxs_name]).astype(np.int32)

res = run_stationary_state_decoder_single_session(
    X_counts,
    state_idxs,
    t_edges=t_edges,
    out_h5_path=snakemake.output[0],
    keys=STATE_KEYS,
    label_names=("tgt_sta", "bgr_sta", "sil_sta"),
    n_splits=5,
    C=1.0,
    max_iter=2000,
    do_balance=True,
    balance_n_per_class=None,   # None -> min count
    do_label_shuffle=True,
    n_label_shuffles=50,
    rng_seed=0,
)