import numpy as np
from sklearn.base import clone
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.metrics import confusion_matrix
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression


def decode_state_from_lfp(X_lfp, y, n_splits=5, random_state=0):
    """
    X_lfp: (N, n_features)
    y:     (N,) integer state labels
    """

    clf = Pipeline([
        ("scaler", StandardScaler()),
        ("logreg", LogisticRegression(
            multi_class="multinomial",
            solver="lbfgs",
            max_iter=1000,
            class_weight="balanced"
        ))
    ])

    # Cross-validated accuracy
    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state)
    acc = cross_val_score(clf, X_lfp, y, cv=cv, scoring="accuracy")

    # Fit once on all data for confusion matrix
    clf.fit(X_lfp, y)
    y_pred = clf.predict(X_lfp)
    cm = confusion_matrix(y, y_pred, labels=np.unique(y))

    return clf, acc, cm, y_pred


def run_decoder_with_shuffles(
    X,
    y,
    n_splits=5,
    n_shuffles=200,
    random_state=0,
):
    """
    Run real decoder + shuffled-label null for one session & one decoder.

    Returns:
        acc_real        : float
        cm_real         : (K, K) array
        y_real          : (N,) array (same as input y)
        y_pred_real     : (N,) array
        acc_shuf        : (n_shuffles,) array
        cm_shuf_mean    : (K, K) array
        cm_shuf_std     : (K, K) array
    """

    # --- REAL DECODER ---
    clf, acc_cv, cm_real, y_pred_real = decode_state_from_lfp(
        X, y, n_splits=n_splits, random_state=random_state
    )
    acc_real = float(np.mean(acc_cv))
    y_real = y.copy()

    # --- SHUFFLES / NULL ---
    n_classes = cm_real.shape[0]
    acc_shuf = np.zeros(n_shuffles, dtype=float)
    cms_shuf = np.zeros((n_shuffles, n_classes, n_classes), dtype=float)

    rng = np.random.default_rng(random_state)
    labels = np.unique(y)

    # define CV once (same strategy for all shuffles)
    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state)

    for s in range(n_shuffles):
        # 1) shuffle labels
        y_shuf = rng.permutation(y)

        # 2) cross-validated accuracy for shuffled labels
        fold_accs = []
        for train_idx, test_idx in cv.split(X, y_shuf):
            clf_fold = clone(clf)
            clf_fold.fit(X[train_idx], y_shuf[train_idx])
            fold_accs.append(clf_fold.score(X[test_idx], y_shuf[test_idx]))
        acc_shuf[s] = np.mean(fold_accs)

        # 3) confusion matrix for shuffled labels (fit on all data)
        clf_full = clone(clf)
        clf_full.fit(X, y_shuf)
        y_pred_shuf = clf_full.predict(X)
        cms_shuf[s] = confusion_matrix(y_shuf, y_pred_shuf, labels=labels)

    cm_shuf_mean = cms_shuf.mean(axis=0)
    cm_shuf_std = cms_shuf.std(axis=0)

    return (
        acc_real,
        cm_real,
        y_real,
        y_pred_real,
        acc_shuf,
        cm_shuf_mean,
        cm_shuf_std,
    )

# --- example usage ---
# acc_real, cm_real, y_real, y_pred_real, acc_shuf, cm_shuf_mean, cm_shuf_std = \
#     run_decoder_with_shuffles(X, y, n_splits=5, n_shuffles=200, random_state=0)

# p-value example:
# p_value = (np.sum(acc_shuf >= acc_real) + 1) / (len(acc_shuf) + 1)
