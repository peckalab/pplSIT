import os, json
import h5py
import numpy as np
from scipy import signal
from scipy.ndimage import median_filter


def head_direction(tl, hd_update_speed=0.04):
    width = 200  # 100 points ~= 1 sec at 100Hz
    kernel = signal.windows.gaussian(width, std=(width) / 7.2)

    x_smooth = np.convolve(tl[:, 1], kernel, "same") / kernel.sum()
    y_smooth = np.convolve(tl[:, 2], kernel, "same") / kernel.sum()

    diff_x = np.diff(x_smooth, axis=0)
    diff_y = np.diff(y_smooth, axis=0)

    hd = -np.arctan2(diff_y, diff_x)
    hd = np.concatenate([np.array([hd[0]]), hd])  # same length as timeline

    # reset idle periods
    idle_idxs = np.where(tl[:, 3] < hd_update_speed)[0]
    if len(idle_idxs) == 0:
        return hd

    crit = np.where(np.diff(idle_idxs) > 1)[0]
    if len(crit) == 0:
        # whole idle stretch
        i1, i2 = idle_idxs[0], idle_idxs[-1]
        hd[i1:i2] = hd[i1 - 1] if i1 > 0 else hd[i1]
        return hd

    idle_periods = []
    idle_periods.append((idle_idxs[0], idle_idxs[crit[0]]))
    for i in range(len(crit) - 1):
        idx_start = idle_idxs[crit[i] + 1]
        idx_end = idle_idxs[crit[i + 1]]
        idle_periods.append((idx_start, idx_end))

    for (i1, i2) in idle_periods:
        hd[i1:i2] = hd[i1 - 1] if i1 > 0 else hd[i1]

    return hd


def build_tgt_matrix(sound_events, trials):
    tgt_start_idxs = []
    tgt_end_idxs = []

    for i in range(len(sound_events) - 1):
        if sound_events[i][1] != 2 and sound_events[i + 1][1] == 2:
            tgt_start_idxs.append(i + 1)
        if sound_events[i][1] == 2 and sound_events[i + 1][1] != 2:
            tgt_end_idxs.append(i)

    if len(tgt_start_idxs) == 0 or len(tgt_end_idxs) == 0:
        return np.zeros((0, 5), dtype=np.int32)

    if tgt_start_idxs[-1] > tgt_end_idxs[-1]:
        tgt_start_idxs = tgt_start_idxs[:-1]
    if len(tgt_start_idxs) == 0:
        return np.zeros((0, 5), dtype=np.int32)

    if tgt_end_idxs[0] < tgt_start_idxs[0]:
        tgt_start_idxs = [0] + tgt_start_idxs

    tgt_start_idxs = np.array(tgt_start_idxs)
    tgt_end_idxs = np.array(tgt_end_idxs)

    tgt_results = np.zeros(len(tgt_start_idxs))
    for idx_tl_success_end in trials[trials[:, 5] == 1][:, 1]:
        idx_succ = np.abs(sound_events[tgt_end_idxs][:, 2] - idx_tl_success_end).argmin()
        tgt_results[idx_succ] = 1

    return np.column_stack(
        [
            tgt_start_idxs,
            tgt_end_idxs,
            sound_events[tgt_start_idxs][:, 2],
            sound_events[tgt_end_idxs][:, 2],
            tgt_results,
        ]
    ).astype(np.int32)


def build_dis_matrix(sound_events, trials, cfg):
    dis_start_idxs = []
    dis_end_idxs = []

    for i in range(len(sound_events) - 1):
        if sound_events[i][1] <= 2 and sound_events[i + 1][1] > 2:
            dis_start_idxs.append(i + 1)
        if sound_events[i][1] > 2 and sound_events[i + 1][1] <= 2:
            dis_end_idxs.append(i)

    if len(dis_start_idxs) == 0 or len(dis_end_idxs) == 0:
        return np.zeros((0, 5), dtype=np.int32)

    if dis_start_idxs[-1] > dis_end_idxs[-1]:
        dis_start_idxs = dis_start_idxs[:-1]
    if len(dis_start_idxs) == 0:
        return np.zeros((0, 5), dtype=np.int32)

    if dis_end_idxs[0] < dis_start_idxs[0]:
        dis_start_idxs = [0] + dis_start_idxs

    dis_start_idxs = np.array(dis_start_idxs)
    dis_end_idxs = np.array(dis_end_idxs)

    dis_results = np.zeros(len(dis_start_idxs))

    fails = trials[trials[:, 5] == 0]
    for idx_tl_fail_start, idx_tl_fail_end in zip(fails[:, 0], fails[:, 1]):
        dis_idxs_before_end = np.where(sound_events[:, 2] <= idx_tl_fail_end)[0]
        if len(dis_idxs_before_end) == 0:
            continue
        last_sound_idx = dis_idxs_before_end[-1]
        if sound_events[last_sound_idx][1] > 2:
            # only if trial was not full duration it is a true fail
            if (idx_tl_fail_end - idx_tl_fail_start) < cfg["experiment"]["trial_duration"] * 100:
                dis_results[np.where(dis_end_idxs == last_sound_idx)[0][0]] = 1
            else:
                # timeout (invalid)
                dis_results[np.where(dis_end_idxs == last_sound_idx)[0][0]] = -1

    return np.column_stack(
        [
            dis_start_idxs,
            dis_end_idxs,
            sound_events[dis_start_idxs][:, 2],
            sound_events[dis_end_idxs][:, 2],
            dis_results,
        ]
    ).astype(np.int32)


def pack_base(
    pos_file,
    ev_file,
    snd_file,
    isl_file,
    cfg_file,
    man_file,
    dst_file,
    drift_coeff=0.000025,
):
    """
    Behavior-only base pack.
    - Always writes raw positions/events/sounds/(islands) and processed datasets.
    - If manual offset is int: apply drift + offset correction to logged sounds (as in original script).
    - If manual offset is dict: DO NOT sync (no ephys inputs here), keep sounds as logged.
      A flag is stored in HDF5 attrs indicating ephys sync was requested.
    """

    with open(cfg_file) as jf:
        parameters = json.load(jf)

    with open(man_file) as jf:
        manual = json.load(jf)

    offset = manual.get("ephys", {}).get("offset", 0)

    # ---------- write HDF5 ----------
    os.makedirs(os.path.dirname(dst_file), exist_ok=True)
    with h5py.File(dst_file, "w") as f:
        raw = f.create_group("raw")
        raw.attrs["parameters"] = json.dumps(parameters)
        raw.attrs["manual"] = json.dumps(manual)

        # -------- save raw CSVs ------------
        ds_names = ["positions", "events", "sounds", "islands"]
        for ds_name, f_path in zip(ds_names, [pos_file, ev_file, snd_file, isl_file]):
            if (not f_path) or (not os.path.exists(f_path)):
                continue
            with open(f_path) as ff:
                headers = ff.readline()
            data = np.loadtxt(f_path, delimiter=",", skiprows=1)
            ds = raw.create_dataset(ds_name, data=data)
            ds.attrs["headers"] = headers

        # read raw data and normalize to session start
        events = np.array(f["raw"]["events"])
        s_start, s_end = events[:, 0][0], events[:, 0][-1]
        events[:, 0] -= s_start

        positions = np.array(f["raw"]["positions"])
        positions[:, 0] -= s_start

        sounds = np.array(f["raw"]["sounds"])
        sounds[:, 0] -= s_start

        # squeeze - if session was interrupted, adjust times to have a continuous timeline
        end_idxs = np.where(events[:, 5] == -1)[0]
        if len(end_idxs) > 1:
            deltas = [events[idx + 1][0] - events[idx][0] for idx in end_idxs[:-1]]

            for df, delta in zip(end_idxs, deltas):
                events[df + 1 :][:, 0] -= delta

            end_idxs_s = np.where(np.diff(sounds[:, 0]) > 20)[0]
            for df, delta in zip(end_idxs_s, deltas):
                sounds[df + 1 :][:, 0] -= delta

            end_idxs_p = np.where(np.diff(positions[:, 0]) > 20)[0]
            for df, delta in zip(end_idxs_p, deltas):
                positions[df + 1 :][:, 0] -= delta

            parameters["experiment"]["timepoints"] = [positions[df + 1][0] for df in end_idxs_p]
            parameters["experiment"]["session_duration"] = positions[-1][0]

        # -------- processed ------------
        proc = f.create_group("processed")
        proc.attrs["parameters"] = json.dumps(parameters)

        # convert timeline to 100 Hz
        time_freq = 100
        s_start, s_end = events[0][0], events[-1][0]
        times = np.arange(s_start, s_end, 1.0 / time_freq)
        pos_at_freq = np.zeros((len(times), 3))

        curr_idx = 0
        for i, t in enumerate(times):
            if (
                curr_idx < len(positions) - 1
                and np.abs(t - positions[:, 0][curr_idx]) > np.abs(t - positions[:, 0][curr_idx + 1])
            ):
                curr_idx += 1
            pos_at_freq[i] = (t, positions[curr_idx][1], positions[curr_idx][2])

        # trials
        t_count = len(np.unique(events[events[:, -1] > 0][:, -2]))
        trials = np.zeros((t_count, 6))
        for i in range(t_count):
            t_start_idx = (np.abs(pos_at_freq[:, 0] - events[2 * i][0])).argmin()
            t_end_idx = (np.abs(pos_at_freq[:, 0] - events[2 * i + 1][0])).argmin()
            state = 1 if events[2 * i + 1][-1] == 1 else 0
            trials[i] = (
                t_start_idx,
                t_end_idx,
                events[2 * i][1],
                events[2 * i][2],
                events[2 * i][3],
                state,
            )

        ds = proc.create_dataset("trial_idxs", data=trials)
        ds.attrs["headers"] = "t_start_idx, t_end_idx, target_x, target_y, target_r, fail_or_success"

        # ---- drift/offset handling (base only) ----
        if isinstance(offset, int):
            drift = s_end * drift_coeff
            sounds[:, 0] = sounds[:, 0] + np.arange(len(sounds)) * drift / len(sounds) + offset / 1000.0
            proc.attrs["sound_time_correction"] = json.dumps(
                {"mode": "manual_drift", "offset_ms": offset, "drift_coeff": drift_coeff}
            )
        else:
            # ephys sync requested, but base pack does not do it
            proc.attrs["sound_time_correction"] = json.dumps(
                {"mode": "no_sync_in_base", "requested": True, "offset_spec": offset}
            )
            raw.attrs["sync_requested"] = True
            raw.attrs["sync_type"] = str(getattr(offset, "get", lambda *_: "unknown")("type", "unknown"))

        # sound_events: map sounds to timeline indices
        sound_events = np.zeros((len(sounds), 3))
        left_idx = 0
        delta = 1e5
        for i in range(len(sounds)):
            while left_idx < len(pos_at_freq) and np.abs(sounds[i][0] - pos_at_freq[:, 0][left_idx]) < delta:
                delta = np.abs(sounds[i][0] - pos_at_freq[:, 0][left_idx])
                left_idx += 1

            sound_events[i] = (sounds[i][0], sounds[i][1], left_idx)
            delta = 1e5

        # remove events outside timeline
        to_keep_idxs = np.where(sound_events[:, 2] < len(times))[0]
        sound_events = sound_events[to_keep_idxs]

        ds = proc.create_dataset("sound_events", data=sound_events)
        ds.attrs["headers"] = "sound_time, sound_id, timeline_idx"

        # smooth + speed + HD
        x_mf = median_filter(pos_at_freq[:, 1], size=200)
        y_mf = median_filter(pos_at_freq[:, 2], size=200)

        width = 100
        kernel = signal.windows.gaussian(width, std=(width) / 7.2)

        x_smooth = np.convolve(x_mf, kernel, "same") / kernel.sum()
        y_smooth = np.convolve(y_mf, kernel, "same") / kernel.sum()

        dx = np.sqrt(np.square(np.diff(x_smooth)) + np.square(np.diff(y_smooth)))
        dt = np.diff(pos_at_freq[:, 0])
        speed = np.concatenate([dx / dt, [dx[-1] / dt[-1]]])

        temp_tl = np.column_stack([pos_at_freq[:, 0], x_smooth, y_smooth, speed])
        hd = head_direction(temp_tl)

        # trial numbers in timeline
        trials_data = np.zeros(len(temp_tl))
        for i, trial in enumerate(trials):
            idx1, idx2 = int(trial[0]), int(trial[1])
            trials_data[idx1:idx2] = i + 1

        # sounds played in timeline
        sound_tl = np.zeros(len(temp_tl))
        curr_sound_idx = 0
        for i in range(len(temp_tl)):
            if curr_sound_idx + 1 >= len(sounds):
                break
            if temp_tl[i][0] > sounds[curr_sound_idx][0]:
                curr_sound_idx += 1
            sound_tl[i] = sounds[curr_sound_idx][1]

        timeline = proc.create_dataset(
            "timeline",
            data=np.column_stack(
                [pos_at_freq[:, 0], x_smooth, y_smooth, speed, hd, trials_data, sound_tl, x_smooth, y_smooth]
            ),
        )
        timeline.attrs["headers"] = "time, x, y, speed, hd, trial_no, sound_ids, x_raw, y_raw"

        # target matrix
        tgt_matrix = build_tgt_matrix(sound_events, trials)
        ds = proc.create_dataset("target_matrix", data=tgt_matrix)
        ds.attrs["headers"] = "sound_idx_start, sound_idx_end, tl_idx_start, tl_idx_end, result"

        # distractor matrix
        if parameters["experiment"].get("distractor_fail", False):
            dis_matrix = build_dis_matrix(sound_events, trials, parameters)
            ds = proc.create_dataset("distractor_matrix", data=dis_matrix)
            ds.attrs["headers"] = "sound_idx_start, sound_idx_end, tl_idx_start, tl_idx_end, result"


# ---- Snakemake entry point ----
pack_base(
    snakemake.input.positions,
    snakemake.input.events,
    snakemake.input.sounds,
    snakemake.input.islands,
    snakemake.input.cfg,
    snakemake.input.manual,
    snakemake.output.base,
    drift_coeff=snakemake.config["pack"]["drift_coeff"],
)
