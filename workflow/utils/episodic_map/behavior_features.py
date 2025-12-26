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
    arena_radius_m=None,
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

    # ---- Derived geometry covariates (center-referenced) ----
    # Distance to arena center (0,0)
    Lc = min(len(neck_x_s), len(neck_y_s), len(t))
    neck_x_s = neck_x_s[:Lc]
    neck_y_s = neck_y_s[:Lc]
    t = t[:Lc]
    pos_valid = pos_valid[:Lc]

    dist_to_center_s = np.sqrt(neck_x_s**2 + neck_y_s**2)

    # Angle from animal position pointing TO the center
    angle_to_center_s = np.arctan2(-neck_y_s, -neck_x_s)

    # HD relative to the center direction (egocentric bearing to center)
    # hd_rel_center = 0 means facing the center
    Lh = min(len(hd_angle_s), len(angle_to_center_s), len(hd_valid), len(pos_valid))
    hd_angle_s = hd_angle_s[:Lh]
    angle_to_center_s = angle_to_center_s[:Lh]
    hd_valid = hd_valid[:Lh]
    pos_valid = pos_valid[:Lh]
    dist_to_center_s = dist_to_center_s[:Lh]
    neck_x_s = neck_x_s[:Lh]
    neck_y_s = neck_y_s[:Lh]
    t = t[:Lh]

    hd_rel_center_s = wrap_to_pi(hd_angle_s - angle_to_center_s)

    # Optional: distance to boundary if arena radius known (circular arena)
    dist_to_boundary_s = None
    if arena_radius_m is not None:
        dist_to_boundary_s = float(arena_radius_m) - dist_to_center_s
        # negative values can happen if tracking drifts outside; clip if you want
        # dist_to_boundary_s = np.clip(dist_to_boundary_s, 0.0, None)

    # ---- Episode features ----
    N_ep = targets_periods_times.shape[0]
    hd_mean = np.full(N_ep, np.nan)
    hd_var = np.full(N_ep, np.nan)
    stillness = np.full(N_ep, np.nan)
    head_turn_rate = np.full(N_ep, np.nan)
    valid_frac = np.full(N_ep, np.nan)
    dist_center_ep = np.full(N_ep, np.nan)              # NEW
    hd_rel_center_mean = np.full(N_ep, np.nan)          # NEW
    hd_rel_center_var = np.full(N_ep, np.nan)           # NEW
    hd_rel_center_valid_frac = np.full(N_ep, np.nan)    # NEW

    dist_boundary_ep = None
    if arena_radius_m is not None:
        dist_boundary_ep = np.full(N_ep, np.nan)

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

        # ---- Center distance (position-only) ----
        pos_ok = np.isfinite(dist_to_center_s) & pos_valid
        if np.any(mask & pos_ok):
            # median is robust to occasional jumps
            dist_center_ep[i] = float(np.nanmedian(dist_to_center_s[mask & pos_ok]))
            if arena_radius_m is not None and dist_to_boundary_s is not None:
                dist_boundary_ep[i] = float(np.nanmedian(dist_to_boundary_s[mask & pos_ok]))

        # ---- HD relative to center (requires both HD and position) ----
        rel_ok = np.isfinite(hd_rel_center_s) & hd_valid & pos_valid
        n_mask = int(np.sum(mask))
        if n_mask > 0:
            hd_rel_center_valid_frac[i] = float(np.sum(mask & rel_ok) / n_mask)

        if np.sum(mask & rel_ok) >= 5:
            th = hd_rel_center_s[mask & rel_ok]
            hd_rel_center_mean[i] = float(circular_mean(th))
            hd_rel_center_var[i]  = float(circular_variance(th))

    session_ts = dict(
        time_s=t,
        hd_angle_rad=hd_angle_s,
        head_ang_vel_rads=head_ang_vel_s,
        body_speed_rms_mps=body_speed_rms_s,
        hd_valid_mask=hd_valid.astype(np.uint8),
        pos_x_m=neck_x_s,
        pos_y_m=neck_y_s,
        pos_valid_mask=pos_valid.astype(np.uint8),
        dist_to_center_m=dist_to_center_s,
        angle_to_center_rad=angle_to_center_s,
        hd_rel_center_rad=hd_rel_center_s,
    )
    if dist_to_boundary_s is not None:
        session_ts["dist_to_boundary_m"] = dist_to_boundary_s

    episode_feats = dict(
        hd_mean_rad=hd_mean,
        hd_var=hd_var,
        stillness_mps=stillness,
        head_turn_rate_rads=head_turn_rate,
        hd_valid_frac=valid_frac,
        targets_periods_times=targets_periods_times,
        dist_to_center_m=dist_center_ep,
        hd_rel_center_mean_rad=hd_rel_center_mean,
        hd_rel_center_var=hd_rel_center_var,
        hd_rel_center_valid_frac=hd_rel_center_valid_frac,
    )
    if dist_boundary_ep is not None:
        episode_feats["dist_to_boundary_m"] = dist_boundary_ep

    return session_ts, episode_feats


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
