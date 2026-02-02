from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import h5py
import numpy as np
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge


# ----------------------------
# Utilities
# ----------------------------

import numpy as np
from typing import Optional

def _build_stationary_windows(
    mask_allowed: np.ndarray,
    win_bins: int,
    hop_bins: Optional[int] = None,
) -> np.ndarray:
    """
    Return window [start, end] indices for contiguous windows length win_bins fully inside mask_allowed.

    IMPORTANT:
      If we return *all* valid starts, we get huge counts because windows overlap (shift by 1 bin).
      Here we downsample starts within each contiguous stationary segment using hop_bins.

    Returns:
      windows: (n_windows, 2) int64 array, each row is [start_idx, end_idx] inclusive.
    """
    if hop_bins is None:
        hop_bins = int(win_bins)
    hop_bins = max(1, int(hop_bins))

    mv = np.asarray(mask_allowed, dtype=np.int8)
    if mv.ndim != 1:
        raise ValueError("mask_allowed must be 1D")
    if win_bins <= 0:
        raise ValueError("win_bins must be > 0")

    # Find contiguous True segments in mask_allowed
    d = np.diff(np.r_[0, mv, 0])
    run_starts = np.where(d == 1)[0]
    run_ends = np.where(d == -1)[0] - 1  # inclusive

    windows = []
    for a, b in zip(run_starts, run_ends):
        run_len = b - a + 1
        if run_len < win_bins:
            continue

        last_start = b - win_bins + 1
        for s in range(a, last_start + 1, hop_bins):
            e = s + win_bins - 1
            windows.append((s, e))

    if len(windows) == 0:
        return np.zeros((0, 2), dtype=np.int64)

    return np.asarray(windows, dtype=np.int64)


def bin_timeseries_to_latent_bins(
    ts_time_s: np.ndarray,
    ts_values: np.ndarray,
    n_bins: int,
    bin_size_s: float,
    method: str = "median",
):
    """
    Aggregate a high-rate time series (e.g. 100 Hz DLC) into latent bins (e.g. 50 ms).

    Parameters
    ----------
    ts_time_s : (N,) float
        Monotonic timestamps in seconds.
    ts_values : (N,) float
        Time series values (may contain NaNs).
    n_bins : int
        Number of latent bins (len(Z_all)).
    bin_size_s : float
        Latent bin size (e.g. 0.05).
    method : str
        'median' or 'mean'.

    Returns
    -------
    binned : (n_bins,) float
        Aggregated values per latent bin (NaN if no samples).
    """
    ts_time_s = np.asarray(ts_time_s).astype(float)
    ts_values = np.asarray(ts_values).astype(float)

    t0 = ts_time_s[0]
    edges = t0 + np.arange(n_bins + 1) * float(bin_size_s)

    binned = np.full(n_bins, np.nan, dtype=float)

    # indices per bin using searchsorted (fast, robust)
    idx0 = np.searchsorted(ts_time_s, edges[:-1], side="left")
    idx1 = np.searchsorted(ts_time_s, edges[1:], side="left")

    for i in range(n_bins):
        a, b = int(idx0[i]), int(idx1[i])
        if b <= a:
            continue
        seg = ts_values[a:b]
        if seg.size == 0:
            continue
        if method == "median":
            binned[i] = np.nanmedian(seg)
        elif method == "mean":
            binned[i] = np.nanmean(seg)
        else:
            raise ValueError(f"Unknown method={method}")

    return binned

def select_window_starts_strided(
    valid_starts: np.ndarray,
    win_bins: int,
    stride_bins: int | None = None,
) -> np.ndarray:
    """
    Convert a dense list of all-valid sliding-window starts into a sparse list
    of starts, enforcing a minimum separation to avoid massive overlap.

    If stride_bins is None -> default stride = win_bins (non-overlapping).
    """
    if valid_starts.size == 0:
        return valid_starts

    if stride_bins is None:
        stride_bins = int(win_bins)

    valid_starts = np.asarray(valid_starts, dtype=np.int64)
    valid_starts.sort()

    chosen = []
    next_allowed = -10**18
    for s in valid_starts:
        if s >= next_allowed:
            chosen.append(int(s))
            next_allowed = int(s) + int(stride_bins)

    return np.asarray(chosen, dtype=np.int64)


def build_stationary_episode_trajectories(
    core_h5: str,
    beh_h5: str,
    *,
    representation: str = "raw",  # raw | resid | resid_stim | resid_ctx_stim
    speed_thresh_mps: float = 0.04,
    episode_len_s: float = 3.0,
    margin_s: float = 0.5,
):
    """
    Adds non-target stationary episode trajectories to existing core HDF5.
    """

    with h5py.File(core_h5, "r+") as f:
        meta = json.loads(f.attrs["meta_json"])
        bin_size_s = float(meta["bin_size_s"])
        win_bins = int(round(episode_len_s / bin_size_s))
        margin_bins = int(round(margin_s / bin_size_s))

        # ---- load latent
        if representation == "raw":
            Z_all = f["latent/Z_all"][...]
        elif representation == "resid":
            Z_all = f["latent/Z_all_resid"][...]
        elif representation == "resid_stim":
            Z_all = f["latent/Z_all_resid_stim"][...]
        elif representation == "resid_ctx_stim":
            Z_all = f["latent/Z_all_resid_ctx_stim"][...]
        else:
            raise ValueError(f"Unknown representation {representation}")

        T, D = Z_all.shape

        # ---- load masks
        mask_tgt = f["states/mask_tgt"][...].astype(bool)
        mask_sta = f["states/mask_sta"][...].astype(bool)

        # ---- build exclusion mask (target ± margin)
        mask_excl = mask_tgt.copy()
        idx = np.where(mask_tgt)[0]
        for i in idx:
            a = max(0, i - margin_bins)
            b = min(T, i + margin_bins + 1)
            mask_excl[a:b] = True

        # ---- load behavior speed
        with h5py.File(beh_h5, "r") as fb:
            sp_hi = fb["dlc_features/session_ts/body_speed_rms_mps"][...].astype(float)
            t_hi  = fb["dlc_features/session_ts/time_s"][...].astype(float)

        # Bin DLC speed to latent bins
        sp = bin_timeseries_to_latent_bins(
            ts_time_s=t_hi,
            ts_values=sp_hi,
            n_bins=T,
            bin_size_s=bin_size_s,
            method="median",
        )

        mask_slow = np.isfinite(sp) & (sp < speed_thresh_mps)
        mask_sta = mask_sta & mask_slow

        # ---- helper to extract & store
        def _extract_and_store(group_name, mask_group):
            allowed = mask_group & (~mask_excl)
            windows = _build_stationary_windows(allowed, win_bins=win_bins, hop_bins=win_bins)
            #windows = _build_stationary_windows(
            #    mask_group, mask_excl, win_bins
            #)
            if len(windows) == 0:
                return

            traj_stack = []
            traj_ptr = []
            row = 0

            for bs, be in windows:
                seg = Z_all[bs:be+1]
                if seg.shape[0] != win_bins:
                    continue
                traj_stack.append(seg)
                traj_ptr.append([row, row + win_bins - 1])
                row += win_bins

            if len(traj_stack) == 0:
                return

            traj_stack = np.vstack(traj_stack).astype(np.float32)
            traj_ptr = np.asarray(traj_ptr, dtype=np.int64)

            g = f.require_group(f"episodes_sta/{group_name}/{representation}")
            g.create_dataset("traj_stack", data=traj_stack, compression="gzip")
            g.create_dataset("traj_ptr", data=traj_ptr, compression="gzip")
            g.create_dataset("bin_windows", data=windows, compression="gzip")

            g.attrs["meta_json"] = json.dumps(dict(
                episode_len_s=episode_len_s,
                episode_win_bins=win_bins,
                speed_thresh_mps=speed_thresh_mps,
                margin_s=margin_s,
                representation=representation,
                n_ep=len(traj_ptr),
            ))

        # ---- ALL STA
        _extract_and_store("all", mask_sta)

        # ---- BGR STA
        if "states/bgr_sta_periods_bin" in f:
            m = periods_to_mask(
                f["states/bgr_sta_periods_bin"][...], T
            )
            _extract_and_store("bgr", mask_sta & m)

        # ---- SIL STA
        if "states/sil_sta_periods_bin" in f:
            m = periods_to_mask(
                f["states/sil_sta_periods_bin"][...], T
            )
            _extract_and_store("sil", mask_sta & m)


def clip_periods(periods_bin: np.ndarray, t_bins: int) -> np.ndarray:
    """Clip periods to [0, t_bins-1] and drop invalid (end < start)."""
    p = np.asarray(periods_bin, dtype=np.int64).copy()
    if p.size == 0:
        return p.reshape(0, 2)
    p[:, 0] = np.clip(p[:, 0], 0, t_bins - 1)
    p[:, 1] = np.clip(p[:, 1], 0, t_bins - 1)
    keep = p[:, 1] >= p[:, 0]
    return p[keep]


def periods_to_mask(periods_bin: np.ndarray, t_bins: int) -> np.ndarray:
    """Create boolean mask length t_bins that is True for bins within any period."""
    mask = np.zeros(t_bins, dtype=bool)
    for s, e in np.asarray(periods_bin, dtype=np.int64):
        mask[s:e + 1] = True
    return mask
