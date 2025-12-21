import numpy as np
import h5py
import json


# --------------------------
# Circular helpers
# --------------------------
def wrap_to_pi(theta):
    return (theta + np.pi) % (2 * np.pi) - np.pi


def circular_mean(theta, w=None):
    theta = np.asarray(theta)
    if w is None:
        w = np.ones_like(theta)
    mask = np.isfinite(theta) & np.isfinite(w) & (w > 0)
    if mask.sum() == 0:
        return np.nan
    s = np.sum(w[mask] * np.sin(theta[mask]))
    c = np.sum(w[mask] * np.cos(theta[mask]))
    return np.arctan2(s, c)


def circular_variance(theta, w=None):
    theta = np.asarray(theta)
    if w is None:
        w = np.ones_like(theta)
    mask = np.isfinite(theta) & np.isfinite(w) & (w > 0)
    if mask.sum() == 0:
        return np.nan
    s = np.sum(w[mask] * np.sin(theta[mask]))
    c = np.sum(w[mask] * np.cos(theta[mask]))
    R = np.sqrt(s**2 + c**2) / np.sum(w[mask])
    return 1.0 - R


# --------------------------
# Smoothing helpers
# --------------------------
def moving_average_nan(x, win):
    x = np.asarray(x, dtype=float)
    if win <= 1:
        return x.copy()
    win = int(win)
    pad = win // 2
    valid = np.isfinite(x).astype(float)
    x0 = np.where(np.isfinite(x), x, 0.0)
    xp = np.pad(x0, (pad, pad), mode="edge")
    vp = np.pad(valid, (pad, pad), mode="edge")
    k = np.ones(win)
    num = np.convolve(xp, k, mode="valid")
    den = np.convolve(vp, k, mode="valid")
    out = num / np.maximum(den, 1e-9)
    out[den < 1e-6] = np.nan
    return out


# --------------------------
# Column helpers
# --------------------------
def _col(cols, name):
    return cols.index(name) if name in cols else None


def get_kp(dlc, cols, base):
    ix = _col(cols, f"{base}_x")
    iy = _col(cols, f"{base}_y")
    il = _col(cols, f"{base}_likelihood")
    if ix is None or iy is None or il is None:
        return None, None, None
    return dlc[:, ix], dlc[:, iy], dlc[:, il]


# --------------------------
# Main function
# --------------------------
def compute_dlc_behavior_features(
    dlc_mat,
    dlc_columns,
    targets_periods_times,
    smooth_size=20,
    p_thr=0.8,
    body_points=("neck", "lower_spine", "tail_base"),
    use_ears_for_hd=True,
):
    """
    Robust DLC feature extraction with global alignment.
    """

    dlc = np.asarray(dlc_mat, float)
    cols = list(dlc_columns)
    t = dlc[:, 0]
    dt = np.nanmedian(np.diff(t))

    # ---- Head direction ----
    nose_x, nose_y, nose_l = get_kp(dlc, cols, "nose")
    if use_ears_for_hd:
        l_x, l_y, l_l = get_kp(dlc, cols, "left_ear")
        r_x, r_y, r_l = get_kp(dlc, cols, "right_ear")
    else:
        l_x, l_y, l_l = get_kp(dlc, cols, "left_eye")
        r_x, r_y, r_l = get_kp(dlc, cols, "right_eye")

    base_x = 0.5 * (l_x + r_x)
    base_y = 0.5 * (l_y + r_y)
    base_l = np.minimum(l_l, r_l)

    hd_valid = (
        (nose_l >= p_thr)
        & (base_l >= p_thr)
        & np.isfinite(nose_x)
        & np.isfinite(base_x)
    )

    vx = nose_x - base_x
    vy = nose_y - base_y
    hd_angle = np.full_like(t, np.nan)
    hd_angle[hd_valid] = np.arctan2(vy[hd_valid], vx[hd_valid])
    hd_angle = wrap_to_pi(hd_angle)

    hd_cos = moving_average_nan(np.cos(hd_angle), smooth_size)
    hd_sin = moving_average_nan(np.sin(hd_angle), smooth_size)
    norm = np.sqrt(hd_cos**2 + hd_sin**2)
    hd_cos /= np.maximum(norm, 1e-9)
    hd_sin /= np.maximum(norm, 1e-9)
    hd_angle_s = np.arctan2(hd_sin, hd_cos)

    dcos = np.gradient(hd_cos, dt)
    dsin = np.gradient(hd_sin, dt)
    head_ang_vel = np.abs(hd_cos * dsin - hd_sin * dcos)
    head_ang_vel_s = moving_average_nan(head_ang_vel, smooth_size)

    # ---- Body speed RMS ----
    speeds = []
    for bp in body_points:
        x, y, l = get_kp(dlc, cols, bp)
        if x is None:
            continue

        # --- ALIGN raw arrays first (critical) ---
        L0 = min(len(x), len(y), len(l), len(t))
        x = x[:L0]
        y = y[:L0]
        l = l[:L0]

        # smooth positions
        x_s = moving_average_nan(x, smooth_size)
        y_s = moving_average_nan(y, smooth_size)

        # velocity & speed (same length as x_s/y_s)
        vx = np.gradient(x_s, dt)
        vy = np.gradient(y_s, dt)
        sp = np.sqrt(vx**2 + vy**2)

        # --- ALIGN sp with likelihood (paranoia-safe) ---
        L1 = min(len(sp), len(l))
        sp = sp[:L1]
        l  = l[:L1]

        # likelihood gating
        sp[l < p_thr] = np.nan
        speeds.append(sp)

    if len(speeds) == 0:
        raise KeyError(f"None of the body_points {body_points} were found in DLC columns.")

    # Stack after aligning to common length across points
    Ls = min(len(s) for s in speeds)
    speeds = np.stack([s[:Ls] for s in speeds], axis=1)  # (Ls, n_pts)

    body_speed_rms = np.sqrt(np.nanmean(speeds**2, axis=1))
    body_speed_rms_s = moving_average_nan(body_speed_rms, smooth_size)

    # ---- Animal position (neck) ----
    neck_x, neck_y, neck_l = get_kp(dlc, cols, "neck")
    if neck_x is None:
        raise KeyError("Missing neck_x/y/likelihood in DLC columns (needed for position).")

    # align raw neck arrays to time (paranoia-safe)
    L0 = min(len(neck_x), len(neck_y), len(neck_l), len(t))
    neck_x = neck_x[:L0]
    neck_y = neck_y[:L0]
    neck_l = neck_l[:L0]

    # smooth
    neck_x_s = moving_average_nan(neck_x, smooth_size)
    neck_y_s = moving_average_nan(neck_y, smooth_size)

    # likelihood gating
    L0 = min(len(neck_x), len(neck_y), len(neck_l), len(t))
    neck_x = neck_x[:L0]; neck_y = neck_y[:L0]; neck_l = neck_l[:L0]

    neck_x_s = moving_average_nan(neck_x, smooth_size)
    neck_y_s = moving_average_nan(neck_y, smooth_size)

    L1 = min(len(neck_x_s), len(neck_y_s), len(neck_l))
    neck_x_s = neck_x_s[:L1]
    neck_y_s = neck_y_s[:L1]
    neck_l   = neck_l[:L1]

    pos_valid = np.isfinite(neck_x_s) & np.isfinite(neck_y_s) & (neck_l >= p_thr)
    neck_x_s = np.where(pos_valid, neck_x_s, np.nan)
    neck_y_s = np.where(pos_valid, neck_y_s, np.nan)


    # ---- GLOBAL ALIGNMENT (critical) ----
    L = min(
        len(t),
        len(hd_angle_s),
        len(hd_valid),
        len(head_ang_vel_s),
        len(body_speed_rms_s),
        len(neck_x_s),
        len(neck_y_s),
        len(pos_valid),
    )

    t = t[:L]
    hd_angle_s = hd_angle_s[:L]
    hd_valid = hd_valid[:L]
    head_ang_vel_s = head_ang_vel_s[:L]
    body_speed_rms_s = body_speed_rms_s[:L]

    neck_x_s = neck_x_s[:L]
    neck_y_s = neck_y_s[:L]
    pos_valid = pos_valid[:L]

    # ---- Episode features ----
    N_ep = targets_periods_times.shape[0]
    hd_mean = np.full(N_ep, np.nan)
    hd_var = np.full(N_ep, np.nan)
    stillness = np.full(N_ep, np.nan)
    head_turn_rate = np.full(N_ep, np.nan)
    valid_frac = np.full(N_ep, np.nan)

    for i, (t0, t1) in enumerate(targets_periods_times):
        mask = (t >= t0) & (t <= t1)
        if mask.sum() < 5:
            continue

        hd_ok = np.isfinite(hd_angle_s) & hd_valid
        valid_frac[i] = float(np.mean(hd_ok[mask]))

        theta_seg = hd_angle_s[mask]
        w = hd_ok[mask].astype(float)

        hd_mean[i] = circular_mean(theta_seg, w)
        hd_var[i] = circular_variance(theta_seg, w)
        stillness[i] = float(np.nanmedian(body_speed_rms_s[mask]))
        head_turn_rate[i] = float(np.nanmedian(head_ang_vel_s[mask]))

    session_ts = dict(
        time_s=t,
        hd_angle_rad=hd_angle_s,
        head_ang_vel_rads=head_ang_vel_s,
        body_speed_rms_mps=body_speed_rms_s,
        hd_valid_mask=hd_valid.astype(np.uint8),
        pos_x_m=neck_x_s,
        pos_y_m=neck_y_s,
        pos_valid_mask=pos_valid.astype(np.uint8),
    )

    episode_feats = dict(
        hd_mean_rad=hd_mean,
        hd_var=hd_var,
        stillness_mps=stillness,
        head_turn_rate_rads=head_turn_rate,
        hd_valid_frac=valid_frac,
    )

    return session_ts, episode_feats

# import numpy as np
# import h5py
# import json


# # --------------------------
# # Circular helpers
# # --------------------------
# def wrap_to_pi(theta):
#     """Wrap angle to [-pi, pi)."""
#     return (theta + np.pi) % (2 * np.pi) - np.pi


# def circular_mean(theta, w=None):
#     """Circular mean of angles theta (radians)."""
#     theta = np.asarray(theta)
#     if w is None:
#         w = np.ones_like(theta, dtype=float)
#     else:
#         w = np.asarray(w, dtype=float)

#     mask = np.isfinite(theta) & np.isfinite(w) & (w > 0)
#     if mask.sum() == 0:
#         return np.nan
#     s = np.sum(w[mask] * np.sin(theta[mask]))
#     c = np.sum(w[mask] * np.cos(theta[mask]))
#     return np.arctan2(s, c)


# def circular_variance(theta, w=None):
#     """
#     Circular variance in [0,1], where 0 is tightly clustered.
#     V = 1 - R, where R is mean resultant length.
#     """
#     theta = np.asarray(theta)
#     if w is None:
#         w = np.ones_like(theta, dtype=float)
#     else:
#         w = np.asarray(w, dtype=float)

#     mask = np.isfinite(theta) & np.isfinite(w) & (w > 0)
#     if mask.sum() == 0:
#         return np.nan

#     s = np.sum(w[mask] * np.sin(theta[mask]))
#     c = np.sum(w[mask] * np.cos(theta[mask]))
#     R = np.sqrt(s**2 + c**2) / np.sum(w[mask])
#     return 1.0 - R


# # --------------------------
# # Smoothing helpers
# # --------------------------
# def moving_average(x, win):
#     """
#     Moving average with edge-padding. Works for 1D arrays.
#     If win<=1 returns x unchanged.
#     """
#     x = np.asarray(x, dtype=float)
#     if win is None or win <= 1:
#         return x.copy()
#     win = int(win)
#     pad = win // 2
#     xp = np.pad(x, (pad, pad), mode="edge")
#     k = np.ones(win, dtype=float) / win
#     return np.convolve(xp, k, mode="valid")


# def moving_average_nan(x, win):
#     """
#     Moving average that ignores NaNs by normalizing with valid counts.
#     """
#     x = np.asarray(x, dtype=float)
#     if win is None or win <= 1:
#         return x.copy()
#     win = int(win)
#     pad = win // 2

#     valid = np.isfinite(x).astype(float)
#     x0 = np.where(np.isfinite(x), x, 0.0)

#     xp = np.pad(x0, (pad, pad), mode="edge")
#     vp = np.pad(valid, (pad, pad), mode="edge")

#     k = np.ones(win, dtype=float)
#     num = np.convolve(xp, k, mode="valid")
#     den = np.convolve(vp, k, mode="valid")
#     out = num / np.maximum(den, 1e-9)
#     out[den < 1e-6] = np.nan
#     return out


# # --------------------------
# # DLC column parsing
# # --------------------------
# def _find_col(cols, name):
#     try:
#         return cols.index(name)
#     except ValueError:
#         return None


# def extract_keypoint_xy_lik(dlc, cols, base_name):
#     """
#     Returns x,y,lik arrays for a keypoint, or (None,None,None) if missing.
#     """
#     ix = _find_col(cols, f"{base_name}_x")
#     iy = _find_col(cols, f"{base_name}_y")
#     il = _find_col(cols, f"{base_name}_likelihood")
#     if ix is None or iy is None or il is None:
#         return None, None, None
#     return dlc[:, ix].astype(float), dlc[:, iy].astype(float), dlc[:, il].astype(float)


# def midpoint(a, b):
#     return 0.5 * (a + b)


# # --------------------------
# # Main computation
# # --------------------------
# def compute_dlc_behavior_features(
#     dlc_mat,
#     dlc_columns,
#     targets_periods_times,
#     smooth_size=20,
#     p_thr=0.8,
#     dt_expected=0.01,
#     body_points=("neck", "lower_spine", "tail_base"),
#     use_ears_for_hd=True,
# ):
#     """
#     Parameters
#     ----------
#     dlc_mat : (T, n_cols) array
#         First column must be absolute time in seconds (100Hz).
#     dlc_columns : list[str]
#         Column names matching dlc_mat.
#     targets_periods_times : (N_ep, 2) array
#         Start/end absolute times (sec) for each target episode.
#     smooth_size : int
#         Moving average window in samples (default 20 -> 200ms at 100Hz).
#     p_thr : float
#         Likelihood threshold for head direction validity.
#     dt_expected : float
#         Expected timestep in seconds; used for sanity/derivatives.
#     body_points : tuple[str]
#         Keypoints used for RMS body speed.
#     use_ears_for_hd : bool
#         If True use ear midpoint as head base, else use eye midpoint.

#     Returns
#     -------
#     session_ts : dict of arrays length T
#     episode_feats : dict of arrays length N_ep
#     """

#     dlc = np.asarray(dlc_mat, dtype=float)
#     cols = list(dlc_columns)

#     # time
#     if cols[0] != "time":
#         # allow first col to be time even if named differently
#         # but user said first column is time
#         pass
#     t = dlc[:, 0].astype(float)
#     T = len(t)

#     # estimate dt for derivatives
#     dt = np.nanmedian(np.diff(t))
#     if not np.isfinite(dt) or dt <= 0:
#         dt = dt_expected

#     # --- head direction: nose - head_base
#     nose_x, nose_y, nose_l = extract_keypoint_xy_lik(dlc, cols, "nose")
#     if nose_x is None:
#         raise KeyError("Missing nose_x/y/likelihood in DLC columns")

#     if use_ears_for_hd:
#         le_x, le_y, le_l = extract_keypoint_xy_lik(dlc, cols, "left_ear")
#         re_x, re_y, re_l = extract_keypoint_xy_lik(dlc, cols, "right_ear")
#         if le_x is None or re_x is None:
#             raise KeyError("Missing ear keypoints for head direction (left_ear/right_ear).")
#         base_x = midpoint(le_x, re_x)
#         base_y = midpoint(le_y, re_y)
#         base_l = np.minimum(le_l, re_l)
#     else:
#         le_x, le_y, le_l = extract_keypoint_xy_lik(dlc, cols, "left_eye")
#         re_x, re_y, re_l = extract_keypoint_xy_lik(dlc, cols, "right_eye")
#         if le_x is None or re_x is None:
#             raise KeyError("Missing eye keypoints for head direction (left_eye/right_eye).")
#         base_x = midpoint(le_x, re_x)
#         base_y = midpoint(le_y, re_y)
#         base_l = np.minimum(le_l, re_l)

#     # validity mask for head direction
#     hd_valid = (nose_l >= p_thr) & (base_l >= p_thr) & np.isfinite(nose_x) & np.isfinite(nose_y) & np.isfinite(base_x) & np.isfinite(base_y)

#     vx = nose_x - base_x
#     vy = nose_y - base_y
#     hd_angle = np.full(T, np.nan, dtype=float)
#     hd_angle[hd_valid] = np.arctan2(vy[hd_valid], vx[hd_valid])
#     hd_angle = wrap_to_pi(hd_angle)

#     # smooth head direction in sin/cos space
#     hd_cos = np.cos(hd_angle)
#     hd_sin = np.sin(hd_angle)

#     # masked smoothing (ignore NaNs from invalid frames)
#     hd_cos_s = moving_average_nan(hd_cos, smooth_size)
#     hd_sin_s = moving_average_nan(hd_sin, smooth_size)

#     # renormalize unit vector
#     norm = np.sqrt(hd_cos_s**2 + hd_sin_s**2)
#     hd_cos_s = hd_cos_s / np.maximum(norm, 1e-9)
#     hd_sin_s = hd_sin_s / np.maximum(norm, 1e-9)
#     hd_angle_s = np.arctan2(hd_sin_s, hd_cos_s)

#     # angular velocity |dtheta/dt| computed from sin/cos derivatives to avoid wrap artifacts
#     dcos = np.gradient(hd_cos_s, dt)
#     dsin = np.gradient(hd_sin_s, dt)
#     # angular speed ≈ |dθ/dt| where dθ = (cos*dsin - sin*dcos)
#     head_ang_vel = np.abs(hd_cos_s * dsin - hd_sin_s * dcos)

#     # smooth angular velocity a bit (optional; keep same window)
#     head_ang_vel_s = moving_average_nan(head_ang_vel, smooth_size)

#     # --- body speed RMS across selected core points
#     speeds = []
#     valid_counts = []

#     for bp in body_points:
#         x, y, l = extract_keypoint_xy_lik(dlc, cols, bp)
#         if x is None:
#             continue

#         # smooth positions
#         x_s = moving_average_nan(x, smooth_size)
#         y_s = moving_average_nan(y, smooth_size)

#         # velocity
#         vx_bp = np.gradient(x_s, dt)
#         vy_bp = np.gradient(y_s, dt)
#         sp = np.sqrt(vx_bp**2 + vy_bp**2)

#         # --- ALIGN LENGTHS (critical fix) ---
#         L = min(len(sp), len(l))
#         sp_ = sp[:L]
#         l_  = l[:L]

#         # likelihood gating
#         good = np.isfinite(sp_) & (l_ >= p_thr)
#         sp_g = np.where(good, sp_, np.nan)

#         speeds.append(sp_g)

#     if len(speeds) == 0:
#         raise KeyError(f"None of the body_points {body_points} were found in DLC columns.")

#     speeds = np.stack(speeds, axis=1)  # (T, n_pts)
#     body_speed_rms = np.sqrt(np.nanmean(speeds**2, axis=1))
#     body_speed_rms_s = moving_average_nan(body_speed_rms, smooth_size)

#     # also (optional) nose speed (often helpful)
#     nose_x_s = moving_average_nan(nose_x, smooth_size)
#     nose_y_s = moving_average_nan(nose_y, smooth_size)
#     nose_vx = np.gradient(nose_x_s, dt)
#     nose_vy = np.gradient(nose_y_s, dt)
#     nose_speed = np.sqrt(nose_vx**2 + nose_vy**2)
#     nose_speed_s = moving_average_nan(nose_speed, smooth_size)

#     # --- Episode-level features
#     t_ep = np.asarray(targets_periods_times, dtype=float)
#     if t_ep.ndim != 2 or t_ep.shape[1] != 2:
#         raise ValueError("targets_periods_times must be (N_ep, 2) start/end times in seconds")

#     N_ep = t_ep.shape[0]
#     hd_mean = np.full(N_ep, np.nan, float)
#     hd_var = np.full(N_ep, np.nan, float)
#     stillness = np.full(N_ep, np.nan, float)
#     head_turn_rate = np.full(N_ep, np.nan, float)
#     valid_frac = np.full(N_ep, np.nan, float)

#     for i in range(N_ep):
#         t0, t1 = t_ep[i]
#         if not (np.isfinite(t0) and np.isfinite(t1) and t1 > t0):
#             continue

#         mask = (t >= t0) & (t <= t1)
#         if mask.sum() < 5:
#             continue

#         # validity for HD within episode
#         # use the original hd_valid plus non-nan after smoothing
#         # --- ALIGN LENGTHS ---
#         L = min(len(hd_angle_s), len(hd_valid), len(mask))
#         hd_angle_s_ = hd_angle_s[:L]
#         hd_valid_   = hd_valid[:L]
#         mask_       = mask[:L]

#         hd_ok = np.isfinite(hd_angle_s_) & hd_valid_

#         vf = float(np.mean(hd_ok[mask_])) if mask_.sum() > 0 else np.nan
#         valid_frac[i] = vf

#         theta_seg = hd_angle_s_[mask_]
#         w = hd_ok[mask_].astype(float)

#         hd_mean[i] = circular_mean(theta_seg, w=w)
#         hd_var[i] = circular_variance(theta_seg, w=w)

#         # stillness: use median of smoothed RMS speed
#         stillness[i] = float(np.nanmedian(body_speed_rms_s[mask]))

#         # head turn rate: median smoothed angular velocity
#         head_turn_rate[i] = float(np.nanmedian(head_ang_vel_s[mask]))

#     session_ts = {
#         "time_s": t,
#         "hd_angle_rad": hd_angle_s,          # smoothed
#         "hd_cos": hd_cos_s,
#         "hd_sin": hd_sin_s,
#         "head_ang_vel_rads": head_ang_vel_s, # smoothed |dθ/dt|
#         "body_speed_rms_mps": body_speed_rms_s,
#         "nose_speed_mps": nose_speed_s,
#         "hd_valid_mask": hd_valid.astype(np.uint8),
#     }

#     episode_feats = {
#         "hd_mean_rad": hd_mean,
#         "hd_var": hd_var,
#         "stillness_mps": stillness,
#         "head_turn_rate_rads": head_turn_rate,
#         "hd_valid_frac": valid_frac,
#     }

#     return session_ts, episode_feats


def save_dlc_features_to_h5(
    h5_path,
    session_ts,
    episode_feats,
    group="dlc_features",
    mode="a",
    overwrite_group=True,
    extra_meta=None,
):
    """
    Save DLC-derived session time series and per-episode features to HDF5.

    Parameters
    ----------
    h5_path : str
        Output HDF5 path (existing or new).
    session_ts : dict[str, np.ndarray]
        Whole-session time series (all arrays length T).
    episode_feats : dict[str, np.ndarray]
        Per-episode features (all arrays length N_ep).
    group : str
        Group name under which to store, default "dlc_features".
    mode : str
        'a' append (default) or 'w' overwrite file.
    overwrite_group : bool
        If True, delete existing group if present.
    extra_meta : dict or None
        Optional metadata dict stored as JSON attribute.
    """
    if extra_meta is None:
        extra_meta = {}

    # Basic checks
    # session_ts arrays should have same length
    ts_lens = [len(np.asarray(v)) for v in session_ts.values() if v is not None]
    if len(ts_lens) == 0 or len(set(ts_lens)) != 1:
        raise ValueError("All session_ts arrays must be non-empty and have the same length.")
    T = ts_lens[0]

    # episode_feats arrays should have same length
    ep_lens = [len(np.asarray(v)) for v in episode_feats.values() if v is not None]
    if len(ep_lens) == 0 or len(set(ep_lens)) != 1:
        raise ValueError("All episode_feats arrays must be non-empty and have the same length.")
    N_ep = ep_lens[0]

    meta = {
        "T_session": int(T),
        "N_episodes": int(N_ep),
        "session_keys": list(session_ts.keys()),
        "episode_keys": list(episode_feats.keys()),
        **extra_meta,
    }

    with h5py.File(h5_path, mode) as h5:
        if group in h5 and overwrite_group:
            del h5[group]
        g = h5.require_group(group)

        # store metadata
        g.attrs["meta_json"] = json.dumps(meta)

        # store session time series
        g_ts = g.require_group("session_ts")
        for k, v in session_ts.items():
            if v is None:
                continue
            arr = np.asarray(v)
            if arr.ndim == 1:
                dtype = np.float64 if arr.dtype.kind == "f" else arr.dtype
                g_ts.create_dataset(k, data=arr.astype(dtype), compression="gzip")
            else:
                g_ts.create_dataset(k, data=arr, compression="gzip")

        # store per-episode features
        g_ep = g.require_group("episode_feats")
        for k, v in episode_feats.items():
            if v is None:
                continue
            arr = np.asarray(v)
            dtype = np.float64 if arr.dtype.kind == "f" else arr.dtype
            g_ep.create_dataset(k, data=arr.astype(dtype), compression="gzip")

    return h5_path
