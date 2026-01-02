# episode_window_generalization.py
import json
import numpy as np
import h5py

from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score


# ----------------------------
# HDF5 loading utilities
# ----------------------------

def load_core_meta(h5_path: str) -> dict:
    with h5py.File(h5_path, "r") as f:
        meta = json.loads(f.attrs.get("meta_json", "{}"))
    return meta


def load_target_episode_tensor_from_core(
    h5_path: str,
    *,
    D: int = 10,
    representation: str = "raw",   # "raw" | "resid"
):
    """
    Returns:
      X_ep: (n_ep, T, D) float32
      meta: dict from meta_json
    """
    if representation not in ("raw", "resid"):
        raise ValueError("representation must be 'raw' or 'resid'")

    with h5py.File(h5_path, "r") as f:
        meta = json.loads(f.attrs.get("meta_json", "{}"))
        T = int(meta.get("episode_win_bins", 0))
        if T <= 0:
            raise ValueError(f"{h5_path}: meta_json missing episode_win_bins")

        traj_ptr = np.asarray(f["episodes/traj_ptr"][...], dtype=np.int64)

        if representation == "raw":
            if "episodes/traj_stack" not in f:
                raise KeyError(f"{h5_path}: missing episodes/traj_stack")
            traj_stack = np.asarray(f["episodes/traj_stack"][...], dtype=np.float32)
        else:
            if "episodes/traj_stack_resid" not in f:
                raise KeyError(f"{h5_path}: missing episodes/traj_stack_resid")
            traj_stack = np.asarray(f["episodes/traj_stack_resid"][...], dtype=np.float32)

    n_ep = traj_ptr.shape[0]
    if n_ep == 0:
        raise RuntimeError(f"{h5_path}: n_ep==0")

    n_pc = traj_stack.shape[1]
    if D > n_pc:
        raise ValueError(f"{h5_path}: requested D={D} but only {n_pc} PCs stored")

    # Build tensor
    X_ep = np.empty((n_ep, T, D), dtype=np.float32)
    for i in range(n_ep):
        s, e = traj_ptr[i]
        seg = traj_stack[s:e+1, :D]
        if seg.shape[0] != T:
            raise ValueError(f"{h5_path}: episode {i} length {seg.shape[0]} != T={T}")
        X_ep[i] = seg

    # Drop episodes that contain NaNs/Infs (should be rare after your trajectories fix)
    ok_ep = np.all(np.isfinite(X_ep), axis=(1, 2))
    if not np.all(ok_ep):
        bad = int(np.sum(~ok_ep))
        print(f"[{h5_path}] WARNING: dropping {bad}/{n_ep} episodes with nonfinite values ({representation})")
        X_ep = X_ep[ok_ep]

    meta["n_ep_loaded"] = int(X_ep.shape[0])
    meta["D_used"] = int(D)
    meta["representation"] = representation
    return X_ep, meta


# ----------------------------
# Windowing + mean subtraction
# ----------------------------

POS_NAMES = ("early", "mid", "late")


def window_start_for_pos(T: int, L: int, pos: str) -> int:
    if pos not in POS_NAMES:
        raise ValueError(f"pos must be one of {POS_NAMES}")
    if L > T:
        raise ValueError(f"L={L} > T={T}")

    if pos == "early":
        return 0
    if pos == "late":
        return T - L
    # mid
    return (T - L) // 2


def bins_for_window(T: int, L: int, pos: str) -> np.ndarray:
    s = window_start_for_pos(T, L, pos)
    return np.arange(s, s + L, dtype=np.int64)


def apply_mean_subtract(
    X_ep: np.ndarray,               # (n_ep, T, D)
    mode: str,                      # "none" | "episode_global" | "train_window"
    train_bins: np.ndarray | None,  # (L,) required if train_window
) -> np.ndarray:
    if mode not in ("none", "episode_global", "train_window"):
        raise ValueError("mean_subtract must be 'none'|'episode_global'|'train_window'")

    if mode == "none":
        return X_ep

    X = np.asarray(X_ep, dtype=np.float32)

    if mode == "episode_global":
        mu = np.mean(X, axis=1, keepdims=True)  # (n_ep,1,D)
        return X - mu

    # train_window
    if train_bins is None:
        raise ValueError("train_bins must be provided for mean_subtract='train_window'")
    mu = np.mean(X[:, train_bins, :], axis=1, keepdims=True)  # (n_ep,1,D)
    return X - mu


# ----------------------------
# Decoder
# ----------------------------

def make_clf(C: float = 1.0, max_iter: int = 4000):
    return Pipeline([
        ("scaler", StandardScaler(with_mean=True, with_std=True)),
        ("clf", LogisticRegression(
            C=C,
            solver="lbfgs",
            l1_ratio=0.0,     # explicitly L2, future-proof
            max_iter=max_iter,
        )),
    ])


def flatten_episode_bins(X_ep: np.ndarray, bins: np.ndarray):
    """
    X_ep: (n_ep, T, D)
    bins: (L,)
    Returns:
      X: (n_ep*L, D)
      y: (n_ep*L,) episode IDs
    """
    n_ep, _, D = X_ep.shape
    L = len(bins)
    X = X_ep[:, bins, :].reshape(n_ep * L, D)
    y = np.repeat(np.arange(n_ep, dtype=np.int64), L)
    return X, y


def decode_trainpos_testpos(
    X_ep: np.ndarray,
    *,
    train_bins: np.ndarray,
    test_bins: np.ndarray,
    C: float = 1.0,
):
    """
    Train on bins in train_bins, test on bins in test_bins.
    Each bin is a sample labeled by episode ID.
    """
    Xtr, ytr = flatten_episode_bins(X_ep, train_bins)
    Xte, yte = flatten_episode_bins(X_ep, test_bins)

    clf = make_clf(C=C)
    clf.fit(Xtr, ytr)
    yhat = clf.predict(Xte)
    acc = float(accuracy_score(yte, yhat))
    return acc


# ----------------------------
# Main per-session runner
# ----------------------------

def run_window_generalization_for_session(
    core_h5_path: str,
    out_h5_path: str,
    *,
    D: int = 10,
    L_s_list=(0.25, 0.5, 1.0, 2.0, 3.0),
    representations=("raw", "resid"),
    mean_subtract_modes=("none", "episode_global", "train_window"),
    C: float = 1.0,
):
    """
    Computes a grid of accuracies for each:
      representation x mean_subtract_mode x window_length_s
    where each entry is a 3x3 matrix over train_pos x test_pos (early/mid/late).

    Saves to out_h5_path.
    """
    meta0 = load_core_meta(core_h5_path)
    bin_size_s = float(meta0.get("bin_size_s", 0.05))
    T = int(meta0.get("episode_win_bins", 0))
    if T <= 0:
        raise ValueError(f"{core_h5_path}: meta_json missing episode_win_bins")

    # Prepare output
    out_meta = {
        "source_core_h5": core_h5_path,
        "bin_size_s": bin_size_s,
        "T_episode_bins": T,
        "D": int(D),
        "L_s_list": [float(x) for x in L_s_list],
        "representations": list(representations),
        "mean_subtract_modes": list(mean_subtract_modes),
        "pos_names": list(POS_NAMES),
        "clf": {"type": "LogisticRegression(multinomial)+StandardScaler", "C": float(C)},
    }

    with h5py.File(out_h5_path, "w") as f_out:
        f_out.attrs["meta_json"] = json.dumps(out_meta)

        g = f_out.create_group("window_generalization")

        # store axis labels
        g.create_dataset("pos_names", data=np.array(POS_NAMES, dtype="S"))
        g.create_dataset("L_s", data=np.asarray(L_s_list, dtype=np.float32))
        g.create_dataset("L_bins", data=np.asarray([int(round(L / bin_size_s)) for L in L_s_list], dtype=np.int64))

        # compute
        for rep in representations:
            X_ep, meta = load_target_episode_tensor_from_core(core_h5_path, D=D, representation=rep)
            n_ep = int(X_ep.shape[0])
            chance = float(1.0 / n_ep) if n_ep > 0 else np.nan

            g_rep = g.create_group(rep)
            g_rep.attrs["n_ep"] = n_ep
            g_rep.attrs["chance"] = chance

            for ms in mean_subtract_modes:
                g_ms = g_rep.create_group(ms)

                for L_s in L_s_list:
                    L_bins = int(round(float(L_s) / bin_size_s))
                    if L_bins < 1:
                        raise ValueError(f"L_s={L_s} gives L_bins={L_bins}")

                    if L_bins > T:
                        # window too long: store NaNs
                        acc_mat = np.full((3, 3), np.nan, dtype=np.float32)
                        starts = np.full(3, -1, dtype=np.int64)
                        ds = g_ms.create_group(f"L_{L_s:g}s")
                        ds.create_dataset("acc", data=acc_mat)
                        ds.create_dataset("starts_bins", data=starts)
                        ds.attrs["L_s"] = float(L_s)
                        ds.attrs["L_bins"] = int(L_bins)
                        ds.attrs["note"] = "L_bins > T"
                        continue

                    # compute starts for early/mid/late
                    starts = np.array([window_start_for_pos(T, L_bins, p) for p in POS_NAMES], dtype=np.int64)

                    # 3x3 matrix train_pos x test_pos
                    acc_mat = np.zeros((3, 3), dtype=np.float32)

                    for i_tr, pos_tr in enumerate(POS_NAMES):
                        tr_bins = bins_for_window(T, L_bins, pos_tr)

                        # apply mean subtraction *per train window* if needed
                        X_use = apply_mean_subtract(
                            X_ep,
                            mode=ms,
                            train_bins=tr_bins if ms == "train_window" else None
                        )

                        for j_te, pos_te in enumerate(POS_NAMES):
                            te_bins = bins_for_window(T, L_bins, pos_te)
                            acc = decode_trainpos_testpos(X_use, train_bins=tr_bins, test_bins=te_bins, C=C)
                            acc_mat[i_tr, j_te] = acc

                    ds = g_ms.create_group(f"L_{L_s:g}s")
                    ds.create_dataset("acc", data=acc_mat, compression="gzip")
                    ds.create_dataset("starts_bins", data=starts, compression="gzip")
                    ds.attrs["L_s"] = float(L_s)
                    ds.attrs["L_bins"] = int(L_bins)
                    ds.attrs["n_ep"] = n_ep
                    ds.attrs["chance"] = chance

        return out_h5_path
