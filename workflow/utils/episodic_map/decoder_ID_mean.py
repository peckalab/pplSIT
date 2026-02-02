import json
import numpy as np
import h5py
from dataclasses import dataclass, asdict
from typing import Dict, Tuple, List, Optional


# -----------------------------
# Utilities
# -----------------------------
def _safe_get_dataset(f: h5py.File, path_candidates: List[str]) -> np.ndarray:
    for p in path_candidates:
        if p in f:
            return f[p][...]
    raise KeyError(f"None of these datasets found: {path_candidates}")


def _load_core_episode_tensor(core_h5: str, subgroup: str = None, representation: str = "raw") -> Tuple[np.ndarray, Dict]:
    """
    Returns:
        X_ep: (n_ep, T, D) float32  -- episode tensor in PCA space
        meta: dict                 -- parsed meta_json
    """
    with h5py.File(core_h5, "r") as hfile:
        if subgroup is not None:
            f = hfile[subgroup]
        else:
            f = hfile
        meta = json.loads(f.attrs.get("meta_json", "{}"))
        win_bins = int(meta.get("episode_win_bins", 0))
        if win_bins <= 0:
            raise ValueError("meta_json missing episode_win_bins or invalid")

        traj_ptr = _safe_get_dataset(f, ["episodes/traj_ptr"]).astype(np.int64)

        # Backward/alias support
        rep_alias = {
            "resid_ctx": "resid",
        }
        representation = rep_alias.get(representation, representation)

        rep_to_candidates = {
            "raw": ["episodes/traj_stack"],
            "resid": ["episodes/traj_stack_resid", "episodes/traj_stack_residual", "episodes/traj_stack_res"],
            "resid_stim": ["episodes/traj_stack_resid_stim"],
            "resid_ctx_stim": ["episodes/traj_stack_resid_ctx_stim"],
        }

        if representation not in rep_to_candidates:
            raise ValueError(
                "representation must be one of "
                f"{tuple(rep_to_candidates.keys())}, got {representation!r}"
            )

        traj_stack = _safe_get_dataset(f, rep_to_candidates[representation]).astype(np.float32)

        # if representation == "raw":
        #     traj_stack = _safe_get_dataset(f, ["episodes/traj_stack"]).astype(np.float32)
        # elif representation == "resid":
        #     # your integrated residualization stores these (per your earlier work)
        #     traj_stack = _safe_get_dataset(
        #         f, ["episodes/traj_stack_resid", "episodes/traj_stack_residual", "episodes/traj_stack_res"]
        #     ).astype(np.float32)
        # else:
        #     raise ValueError("representation must be 'raw' or 'resid'")

    n_ep = traj_ptr.shape[0]
    D = traj_stack.shape[1]
    X_ep = np.empty((n_ep, win_bins, D), dtype=np.float32)

    for i in range(n_ep):
        s, e = traj_ptr[i]
        seg = traj_stack[s:e+1]
        if seg.shape[0] != win_bins:
            raise ValueError(f"Episode {i} length mismatch: got {seg.shape[0]}, expected {win_bins}")
        X_ep[i] = seg

    return X_ep, meta


def _make_window_indices(
    T: int,
    bin_size_s: float,
    mean_window: str,
    *,
    strip_first_s: float = 0.0,
    early_s: float = 2.0,
    late_s: float = 2.0,
) -> np.ndarray:
    """
    Returns indices (1D array) into [0..T-1] used for computing episode mean.
    """
    strip = int(round(strip_first_s / bin_size_s))
    strip = max(0, min(strip, T-1))

    if mean_window == "whole":
        idx = np.arange(strip, T, dtype=np.int64)

    elif mean_window == "early":
        n = int(round(early_s / bin_size_s))
        n = max(1, min(n, T - strip))
        idx = np.arange(strip, strip + n, dtype=np.int64)

    elif mean_window == "late":
        n = int(round(late_s / bin_size_s))
        n = max(1, min(n, T - strip))
        idx = np.arange(T - n, T, dtype=np.int64)

    else:
        raise ValueError("mean_window must be one of: 'whole', 'early', 'late'")

    # Ensure at least 2 bins so we can split even/odd into two views
    if idx.size < 2:
        raise ValueError(f"mean_window={mean_window} produced <2 bins after strip. Adjust params.")
    return idx


def _split_even_odd(idx: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Split indices into two disjoint views: even positions vs odd positions in idx.
    """
    # even/odd positions within the idx array (not absolute bin parity)
    a = idx[0::2]
    b = idx[1::2]
    if a.size == 0 or b.size == 0:
        # fallback: if one side empty, make a half/half split
        mid = idx.size // 2
        a = idx[:mid]
        b = idx[mid:]
    if a.size == 0 or b.size == 0:
        raise ValueError("Unable to split bins into two non-empty views.")
    return a, b


def _zscore_train_apply(train: np.ndarray, test: np.ndarray, eps: float = 1e-12) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Z-score features columnwise using train stats only.
    """
    mu = np.nanmean(train, axis=0)
    sd = np.nanstd(train, axis=0)
    sd = np.where(sd < eps, 1.0, sd)
    return (train - mu) / sd, (test - mu) / sd, mu, sd


def _nn_identification_accuracy(trainX: np.ndarray, testX: np.ndarray) -> float:
    """
    1-NN identification: for each episode i, predict label as nearest train prototype.
    Assumes trainX and testX aligned to same episode ordering.
    """
    # pairwise squared Euclidean distances: (n_test, n_train)
    # (x - y)^2 = x^2 + y^2 - 2xy
    x2 = np.sum(testX**2, axis=1, keepdims=True)
    y2 = np.sum(trainX**2, axis=1, keepdims=True).T
    d2 = x2 + y2 - 2.0 * (testX @ trainX.T)
    pred = np.argmin(d2, axis=1)
    true = np.arange(testX.shape[0])
    return float(np.mean(pred == true))


# -----------------------------
# Main API
# -----------------------------
@dataclass
class MeanOnlyDecodeSpec:
    representation: str              # "raw" | "resid" | "resid_stim" | "resid_ctx_stim"
    mean_window: str                # "whole" | "early" | "late"
    strip_first_s: float = 0.0
    early_s: float = 2.0
    late_s: float = 2.0
    zscore_across_episodes: bool = True
    view_split: str = "even_odd"     # currently only even_odd supported


def run_episode_mean_only_decoding_single_session(
    core_h5_path: str,
    out_h5_path: str,
    *,
    D_use: int = 10,
    subgroup: str = None,
    specs: Optional[List[MeanOnlyDecodeSpec]] = None,
) -> str:
    """
    Computes mean-only episode ID decoding for all spec combinations and writes to HDF5.

    Output HDF5 structure:
      /meta_json
      /results/acc               (n_specs,)
      /results/chance            (n_specs,)
      /results/n_ep              (n_specs,)
      /results/representation    (n_specs,) as fixed-length bytes
      /results/mean_window       (n_specs,) as fixed-length bytes
      /results/zscore            (n_specs,) uint8
      /results/strip_first_s     (n_specs,) float32
      /results/early_s           (n_specs,) float32
      /results/late_s            (n_specs,) float32
    """
    if specs is None:
        specs = []
        for rep in ["raw", "resid", "resid_stim", "resid_ctx_stim"]:
            for mw in ["whole", "early", "late"]:
                for z in [False, True]:
                    specs.append(MeanOnlyDecodeSpec(representation=rep, mean_window=mw, zscore_across_episodes=z))

    all_rows = []

    for spec in specs:
        X_ep, meta = _load_core_episode_tensor(core_h5_path, subgroup=subgroup, representation=spec.representation)
        n_ep, T, D = X_ep.shape
        if D_use > D:
            raise ValueError(f"D_use={D_use} > stored PCs={D}")

        bin_size_s = float(meta.get("bin_size_s", 0.05))

        idx = _make_window_indices(
            T, bin_size_s, spec.mean_window,
            strip_first_s=spec.strip_first_s,
            early_s=spec.early_s,
            late_s=spec.late_s,
        )

        if spec.view_split != "even_odd":
            raise ValueError("Only view_split='even_odd' is implemented for now.")
        idx_a, idx_b = _split_even_odd(idx)

        # One prototype per episode, per view
        mu_a = X_ep[:, idx_a, :D_use].mean(axis=1)  # (n_ep, D_use)
        mu_b = X_ep[:, idx_b, :D_use].mean(axis=1)  # (n_ep, D_use)

        # Optional z-score across episodes per PC (fit on train prototypes)
        if spec.zscore_across_episodes:
            mu_a_z, mu_b_z, z_mu, z_sd = _zscore_train_apply(mu_a, mu_b)
            trainX, testX = mu_a_z, mu_b_z
        else:
            trainX, testX = mu_a, mu_b

        acc = _nn_identification_accuracy(trainX, testX)
        chance = 1.0 / float(n_ep)

        all_rows.append({
            "acc": acc,
            "chance": chance,
            "n_ep": int(n_ep),
            "representation": spec.representation,
            "mean_window": spec.mean_window,
            "zscore": bool(spec.zscore_across_episodes),
            "strip_first_s": float(spec.strip_first_s),
            "early_s": float(spec.early_s),
            "late_s": float(spec.late_s),
            "D_use": int(D_use),
        })

    # Write HDF5
    def _as_fixed_str(arr: List[str], L: int = 16) -> np.ndarray:
        return np.array([s.encode("utf-8") for s in arr], dtype=f"|S{L}")

    meta_out = {
        "source_core_h5": core_h5_path,
        "D_use": int(D_use),
        "n_specs": int(len(all_rows)),
        "note": "Episode-mean-only decoding via two-view (even/odd bins) 1-NN identification.",
    }

    with h5py.File(out_h5_path, "w") as f:
        f.attrs["meta_json"] = json.dumps(meta_out)

        g = f.create_group("results")

        g.create_dataset("acc", data=np.array([r["acc"] for r in all_rows], dtype=np.float32))
        g.create_dataset("chance", data=np.array([r["chance"] for r in all_rows], dtype=np.float32))
        g.create_dataset("n_ep", data=np.array([r["n_ep"] for r in all_rows], dtype=np.int32))
        g.create_dataset("D_use", data=np.array([r["D_use"] for r in all_rows], dtype=np.int16))

        g.create_dataset("representation", data=_as_fixed_str([r["representation"] for r in all_rows], L=20))
        g.create_dataset("mean_window", data=_as_fixed_str([r["mean_window"] for r in all_rows], L=8))
        g.create_dataset("zscore", data=np.array([int(r["zscore"]) for r in all_rows], dtype=np.uint8))

        g.create_dataset("strip_first_s", data=np.array([r["strip_first_s"] for r in all_rows], dtype=np.float32))
        g.create_dataset("early_s", data=np.array([r["early_s"] for r in all_rows], dtype=np.float32))
        g.create_dataset("late_s", data=np.array([r["late_s"] for r in all_rows], dtype=np.float32))

    return out_h5_path


# -----------------------------
# Convenience wrapper
# -----------------------------
def run_episode_mean_only_decoding_default(
    core_h5_path: str,
    out_h5_path: str,
    *,
    D_use: int = 10,
    strip_first_s: float = 0.5,
    early_s: float = 2.0,
    late_s: float = 2.0,
) -> str:
    """
    Default grid:
      representation ∈ {raw,resid}
      mean_window ∈ {whole,early,late}
      zscore ∈ {False,True}
    """
    specs = []
    for rep in ["raw", "resid", "resid_stim", "resid_ctx_stim"]:
        for mw in ["whole", "early", "late"]:
            for z in [False, True]:
                specs.append(MeanOnlyDecodeSpec(
                    representation=rep, mean_window=mw,
                    strip_first_s=strip_first_s, early_s=early_s, late_s=late_s,
                    zscore_across_episodes=z
                ))
    return run_episode_mean_only_decoding_single_session(
        core_h5_path, out_h5_path, D_use=D_use, specs=specs
    )
