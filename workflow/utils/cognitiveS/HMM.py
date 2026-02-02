import numpy as np
import h5py
from sklearn.decomposition import PCA
from sklearn.preprocessing import RobustScaler
from hmmlearn.hmm import GaussianHMM
from scipy.ndimage import gaussian_filter1d


# ============================================================
# -------------------- PREPROCESSING -------------------------
# ============================================================

def robust_zscore(X, eps=1e-6):
    """
    Robust z-score per neuron (row-wise).
    X: (neurons, time)
    """
    scaler = RobustScaler(with_centering=True, with_scaling=True)
    Xz = scaler.fit_transform(X.T).T
    return Xz


def smooth_activity(X, sigma_bins=None):
    """
    Optional temporal smoothing (Gaussian).
    X: (neurons, time)
    sigma_bins: None or float (in bins)
    """
    if sigma_bins is None or sigma_bins <= 0:
        return X
    return gaussian_filter1d(X, sigma=sigma_bins, axis=1, mode="nearest")


def subtract_slow_drift(X, sigma_bins=None):
    """
    Subtract very slow drift (running Gaussian mean).
    sigma_bins: None or large (e.g. 600–1200 bins = 30–60 s at 50 ms)
    """
    if sigma_bins is None or sigma_bins <= 0:
        return X
    slow = gaussian_filter1d(X, sigma=sigma_bins, axis=1, mode="nearest")
    return X - slow


# ============================================================
# ------------------ ENSEMBLE FEATURES -----------------------
# ============================================================

def compute_ensemble_scores(Xz, ensemble_labels, min_units=3):
    """
    Compute ensemble score per ensemble using PC1.

    Xz: (neurons, time) z-scored
    ensemble_labels: (neurons,) int labels (0 = ignore)
    Returns:
        features: (time, n_ensembles)
        meta: dict with PCA loadings, unit indices, etc.
    """
    ensembles = np.unique(ensemble_labels)
    ensembles = ensembles[ensembles > 0]

    features = []
    meta = {}

    for g in ensembles:
        idx = np.where(ensemble_labels == g)[0]
        if len(idx) < min_units:
            continue

        Xg = Xz[idx, :].T  # (time, units)

        pca = PCA(n_components=1)
        score = pca.fit_transform(Xg).squeeze()  # (time,)

        features.append(score)
        meta[f"ensemble_{g}"] = {
            "unit_indices": idx,
            "pc1_loadings": pca.components_[0],
            "explained_variance": pca.explained_variance_ratio_[0],
        }

    features = np.column_stack(features)  # (time, n_ensembles)
    return features, meta


def add_global_features(Xz):
    """
    Add global population rate + global PC1.
    """
    global_rate = Xz.mean(axis=0)

    pca = PCA(n_components=1)
    global_pc1 = pca.fit_transform(Xz.T).squeeze()

    return np.column_stack([global_rate, global_pc1]), {
        "global_pc1_variance": pca.explained_variance_ratio_[0]
    }


# ============================================================
# ---------------------- STICKY HMM --------------------------
# ============================================================

def make_sticky_transition_matrix(K, stickiness=10.0):
    """
    Create a diagonal-biased transition matrix prior.
    """
    A = np.ones((K, K))
    A *= 1.0
    np.fill_diagonal(A, stickiness)
    A /= A.sum(axis=1, keepdims=True)
    return A


def fit_sticky_hmm(features, K=4, stickiness=10.0,
                   n_iter=200, random_state=0):
    """
    Fit Gaussian HMM with sticky transition initialization.
    """
    hmm = GaussianHMM(
        n_components=K,
        covariance_type="full",
        n_iter=n_iter,
        random_state=random_state,
        verbose=False,
    )

    # initialize means with k-means-like PCA slicing
    hmm.means_ = np.random.randn(K, features.shape[1])

    hmm.transmat_ = make_sticky_transition_matrix(K, stickiness)
    hmm.startprob_ = np.ones(K) / K

    hmm.fit(features)

    states = hmm.predict(features)
    posteriors = hmm.predict_proba(features)

    return hmm, states, posteriors


# ============================================================
# ---------------------- MAIN PIPELINE -----------------------
# ============================================================

def run_sticky_hmm_pipeline(
    spike_counts,          # (neurons, time)
    ensemble_labels,       # (neurons,)
    bin_size_ms=50,
    smooth_bins=None,      # e.g. 2–20
    drift_bins=None,       # e.g. 600–1200
    K=4,
    stickiness=10.0,
    output_h5="sticky_hmm_results.h5",
):
    """
    Full Tier-A v1 pipeline.
    """

    # ---------- Preprocess ----------
    X = spike_counts.astype(float)

    X = robust_zscore(X)
    X = smooth_activity(X, smooth_bins)
    X = subtract_slow_drift(X, drift_bins)

    # ---------- Ensemble features ----------
    ensemble_feats, ensemble_meta = compute_ensemble_scores(X, ensemble_labels)

    global_feats, global_meta = add_global_features(X)

    features = np.column_stack([ensemble_feats, global_feats])

    # ---------- Fit HMM ----------
    hmm, states, posteriors = fit_sticky_hmm(
        features,
        K=K,
        stickiness=stickiness,
    )

    # ---------- Save ----------
    with h5py.File(output_h5, "w") as f:
        f.create_dataset("features", data=features)
        f.create_dataset("states", data=states)
        f.create_dataset("posteriors", data=posteriors)

        f.create_dataset("transmat", data=hmm.transmat_)
        f.create_dataset("means", data=hmm.means_)
        f.create_dataset("covars", data=hmm.covars_)

        f.attrs["bin_size_ms"] = bin_size_ms
        f.attrs["K"] = K
        f.attrs["stickiness"] = stickiness
        f.attrs["smooth_bins"] = -1 if smooth_bins is None else smooth_bins
        f.attrs["drift_bins"] = -1 if drift_bins is None else drift_bins

        meta_grp = f.create_group("ensemble_meta")
        for k, v in ensemble_meta.items():
            g = meta_grp.create_group(k)
            g.create_dataset("unit_indices", data=v["unit_indices"])
            g.create_dataset("pc1_loadings", data=v["pc1_loadings"])
            g.attrs["explained_variance"] = v["explained_variance"]

        global_grp = f.create_group("global_meta")
        for k, v in global_meta.items():
            global_grp.attrs[k] = v

    return {
        "hmm": hmm,
        "states": states,
        "posteriors": posteriors,
        "features": features,
    }
