import json
from dataclasses import dataclass
from typing import Dict, Optional

import numpy as np
import h5py


# ----------------------------
# Helpers: windows + PCA
# ----------------------------

def _read_bin_size_s(f: h5py.File, default: float = 0.05) -> float:
    if "meta_json" in f.attrs:
        try:
            meta = json.loads(f.attrs["meta_json"])
            if "bin_size_s" in meta:
                return float(meta["bin_size_s"])
        except Exception:
            pass
    if "bin_size_s" in f.attrs:
        return float(f.attrs["bin_size_s"])
    return float(default)


def _get_expected_win_bins(f: h5py.File, key: str) -> Optional[int]:
    if "meta_json" not in f.attrs:
        return None
    try:
        meta = json.loads(f.attrs["meta_json"])
    except Exception:
        return None

    if key == "target":
        for k in ["episode_win_bins", "target_win_bins", "target_episode_win_bins"]:
            if k in meta:
                return int(meta[k])

    if key in ["all_sta", "bgr_sta", "sil_sta"]:
        for k in ["sta_win_bins", "stationary_win_bins", "control_sta_win_bins"]:
            if k in meta:
                return int(meta[k])

    return None


def _load_bin_windows(f: h5py.File, dset_path: str) -> np.ndarray:
    if dset_path not in f:
        raise KeyError(f"Missing dataset: {dset_path}")
    w = f[dset_path][...]
    w = np.asarray(w, dtype=np.int64)
    if w.ndim != 2 or w.shape[1] != 2:
        raise ValueError(f"{dset_path}: expected shape (n_ep,2), got {w.shape}")
    return w


def _detect_inclusive_end(
    windows: np.ndarray,
    n_bins: int,
    expected_win_bins: Optional[int] = None,
    dset_attrs: Optional[h5py.AttributeManager] = None,
) -> bool:
    if dset_attrs is not None:
        for k in ["end_inclusive", "end_is_inclusive", "inclusive_end"]:
            if k in dset_attrs:
                return bool(dset_attrs[k])

    s = windows[:, 0]
    e = windows[:, 1]

    if np.any(e == n_bins):
        return False

    if expected_win_bins is not None and expected_win_bins > 0:
        len_excl = e - s
        len_incl = e - s + 1
        mad_excl = np.median(np.abs(len_excl - expected_win_bins))
        mad_incl = np.median(np.abs(len_incl - expected_win_bins))
        if mad_incl < mad_excl:
            return True
        if mad_excl < mad_incl:
            return False
        return False

    len_excl = e - s
    len_incl = e - s + 1
    bad_excl = np.mean(len_excl <= 0)
    bad_incl = np.mean(len_incl <= 0)
    if bad_incl < bad_excl:
        return True
    return False


def _episode_center_s(starts: np.ndarray, ends: np.ndarray, bin_size_s: float, end_inclusive: bool) -> np.ndarray:
    starts = starts.astype(np.float64)
    ends = ends.astype(np.float64)
    if end_inclusive:
        centers_bin = 0.5 * (starts + ends)
    else:
        centers_bin = 0.5 * (starts + (ends - 1.0))
    return centers_bin * float(bin_size_s)


def _compute_mu_matrix(
    spike_counts: np.ndarray,
    windows: np.ndarray,
    end_inclusive: bool,
) -> np.ndarray:
    sc = np.asarray(spike_counts, dtype=np.float32)
    if sc.ndim != 2:
        raise ValueError(f"spike_counts must be 2D (n_units,n_bins), got {sc.shape}")
    n_units, n_bins = sc.shape

    starts = windows[:, 0].astype(np.int64)
    ends = windows[:, 1].astype(np.int64)

    MU = np.empty((windows.shape[0], n_units), dtype=np.float32)

    for i, (s, e) in enumerate(zip(starts, ends)):
        e2 = (e + 1) if end_inclusive else e
        if s < 0 or e2 > n_bins or e2 <= s:
            raise ValueError(f"Bad window {i}: start={s}, end={e} (inclusive={end_inclusive}), n_bins={n_bins}")
        seg = sc[:, s:e2]
        MU[i] = np.mean(seg, axis=1)

    return MU


def _pca_svd(X: np.ndarray, n_components: int = 10) -> Dict[str, np.ndarray]:
    X = np.asarray(X, dtype=np.float64)
    n_samples, n_features = X.shape
    k = int(min(n_components, n_samples, n_features))
    if k <= 0:
        raise ValueError("n_components too small for data shape")

    U, S, Vt = np.linalg.svd(X, full_matrices=False)
    U = U[:, :k]
    S = S[:k]
    Vt = Vt[:k, :]

    denom = max(1, n_samples - 1)
    explained_variance = (S ** 2) / denom

    # total variance in covariance space
    # (compute from all singular values without recomputing full SVD)
    # Here: use Frobenius norm identity: sum(S_all^2)/(n-1) equals total variance.
    # We already have S_all in SVD; need full singular values:
    S_all = np.linalg.svd(X, compute_uv=False)
    total_var = float(np.sum((S_all ** 2) / denom))
    explained_variance_ratio = explained_variance / total_var if total_var > 0 else np.zeros_like(explained_variance)

    scores = U * S

    return dict(
        components=Vt.astype(np.float32),
        scores=scores.astype(np.float32),
        explained_variance=explained_variance.astype(np.float32),
        explained_variance_ratio=explained_variance_ratio.astype(np.float32),
    )


# ----------------------------
# Main: compute+save drift PCs
# ----------------------------

@dataclass
class DriftPCSpec:
    n_components: int = 10
    zscore_units: bool = True
    overwrite: bool = True


def compute_and_save_drift_pcs_single_session(
    traj_h5_path: str,
    out_h5_path: str,
    spike_counts: np.ndarray,  # (n_units, n_bins)
    episode_sets: Optional[Dict[str, str]] = None,
    out_group_root: str = "drift_pcs",
    spec: DriftPCSpec = DriftPCSpec(),
) -> str:
    """
    Compute drift PCs from per-episode mean firing (MU) for multiple episode sets.

    Reads episode windows from traj_h5_path, writes results to out_h5_path under out_group_root.
    """
    if episode_sets is None:
        episode_sets = {
            "target":   "episodes/target_bin_windows",
            "all_sta":  "episodes_sta/all_sta/raw/bin_windows",
            "bgr_sta":  "episodes_sta/bgr_sta/raw/bin_windows",
            "sil_sta":  "episodes_sta/sil_sta/raw/bin_windows",
        }

    sc = np.asarray(spike_counts)
    if sc.ndim != 2:
        raise ValueError(f"spike_counts must be 2D, got {sc.shape}")
    n_units, n_bins = sc.shape

    with h5py.File(traj_h5_path, "r") as fin, h5py.File(out_h5_path, "a") as fout:
        bin_size_s = _read_bin_size_s(fin, default=0.05)

        # create/clear root in output
        if out_group_root in fout and spec.overwrite:
            del fout[out_group_root]
        root = fout.require_group(out_group_root)

        root.attrs["source_traj_h5"] = str(traj_h5_path)
        root.attrs["bin_size_s"] = float(bin_size_s)
        root.attrs["n_units"] = int(n_units)
        root.attrs["n_bins"] = int(n_bins)
        root.attrs["zscore_units"] = bool(spec.zscore_units)
        root.attrs["n_components_requested"] = int(spec.n_components)

        for set_name, win_path in episode_sets.items():
            w = _load_bin_windows(fin, win_path)

            expected_win = _get_expected_win_bins(fin, key=set_name)
            dset_attrs = fin[win_path].attrs if win_path in fin else None
            end_inclusive = _detect_inclusive_end(w, n_bins=n_bins, expected_win_bins=expected_win, dset_attrs=dset_attrs)

            MU = _compute_mu_matrix(sc, w, end_inclusive=end_inclusive)
            n_ep = MU.shape[0]

            mu_unit = MU.mean(axis=0, keepdims=True)
            X = MU - mu_unit
            if spec.zscore_units:
                sd_unit = MU.std(axis=0, ddof=1, keepdims=True)
                sd_unit = np.where(sd_unit > 0, sd_unit, 1.0)
                X = X / sd_unit
            else:
                sd_unit = np.ones((1, n_units), dtype=np.float32)

            pca = _pca_svd(X, n_components=spec.n_components)

            starts = w[:, 0].astype(np.int64)
            ends = w[:, 1].astype(np.int64)
            centers_s = _episode_center_s(starts, ends, bin_size_s=bin_size_s, end_inclusive=end_inclusive).astype(np.float32)

            g = root.require_group(set_name)
            if spec.overwrite:
                for k in list(g.keys()):
                    del g[k]

            g.attrs["source_windows_path"] = win_path
            g.attrs["end_inclusive"] = bool(end_inclusive)
            if expected_win is not None:
                g.attrs["expected_win_bins"] = int(expected_win)

            g.create_dataset("MU", data=MU.astype(np.float32), compression="gzip")
            g.create_dataset("episode_start_bin", data=starts, compression="gzip")
            g.create_dataset("episode_end_bin", data=ends, compression="gzip")
            g.create_dataset("episode_center_s", data=centers_s, compression="gzip")

            g.create_dataset("unit_mean", data=mu_unit.astype(np.float32).squeeze(0), compression="gzip")
            g.create_dataset("unit_std", data=sd_unit.astype(np.float32).squeeze(0), compression="gzip")

            g_pca = g.require_group("pca")
            g_pca.create_dataset("components", data=pca["components"], compression="gzip")
            g_pca.create_dataset("scores", data=pca["scores"], compression="gzip")
            g_pca.create_dataset("explained_variance", data=pca["explained_variance"], compression="gzip")
            g_pca.create_dataset("explained_variance_ratio", data=pca["explained_variance_ratio"], compression="gzip")

            g_pca.attrs["n_ep"] = int(n_ep)
            g_pca.attrs["n_units"] = int(n_units)
            g_pca.attrs["n_components"] = int(pca["components"].shape[0])

        fout.flush()

    return out_h5_path
