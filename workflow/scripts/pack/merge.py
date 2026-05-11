import os, json, shutil
import numpy as np
import h5py

# --- small utilities copied from pack_base.py (keep consistent) ---

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
    # trials[:,5] == 1 => success; trials[:,1] => trial end index
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
            # find distractor period containing last_sound_idx
            match = np.where((dis_start_idxs <= last_sound_idx) & (last_sound_idx <= dis_end_idxs))[0]

            if len(match) == 0:
                # no matching distractor period found; skip or warn
                continue

            dis_idx = match[-1]

            if (idx_tl_fail_end - idx_tl_fail_start) < cfg["experiment"]["trial_duration"] * 100:
                dis_results[dis_idx] = 1
            else:
                # timeout
                dis_results[dis_idx] = -1

    return np.column_stack(
        [
            dis_start_idxs,
            dis_end_idxs,
            sound_events[dis_start_idxs][:, 2],
            sound_events[dis_end_idxs][:, 2],
            dis_results,
        ]
    ).astype(np.int32)

def nearest_tl_index(times, sound_t):
    # fast nearest index via searchsorted
    # times is sorted timeline time vector at 100 Hz
    i = np.searchsorted(times, sound_t)
    if i <= 0:
        return 0
    if i >= len(times):
        return len(times) - 1
    return i if (sound_t - times[i-1]) >= (times[i] - sound_t) else (i - 1)


# ------------------ main ------------------

base_path = snakemake.input["base"]
sync_none_path = snakemake.input["sync_none"]
# sync_ephys may be absent, [] (Snakemake optional), or a string
sync_ephys_in = snakemake.input.get("sync_ephys", [])
manual_path = snakemake.params.get("manual", None)
out_path  = snakemake.output["meta"]

def _as_one_path(x):
    if x is None:
        return None
    if isinstance(x, str):
        return x
    # snakemake can pass lists for optional inputs
    if isinstance(x, (list, tuple)):
        return x[0] if len(x) else None
    return None

sync_ephys_path = _as_one_path(sync_ephys_in)

os.makedirs(os.path.dirname(out_path), exist_ok=True)
shutil.copy2(base_path, out_path)

# Choose which sync file to use (runtime decision; safe)
use_ephys = False
requested = False
requested_type = "unknown"

if manual_path is not None and os.path.exists(manual_path):
    with open(manual_path, "r") as f:
        manual = json.load(f)
    off = manual.get("ephys", {}).get("offset", 0)
    requested = isinstance(off, dict)
    if requested:
        requested_type = off.get("type", "unknown")

# Prefer ephys sync if requested and available
if requested and sync_ephys_path and os.path.exists(sync_ephys_path):
    sync_path = sync_ephys_path
    use_ephys = True
else:
    sync_path = sync_none_path

with h5py.File(sync_path, "r") as sync:
    mode = sync.attrs.get("mode", "no_sync")
    ev_detected = sync["events_detected"][()] if "events_detected" in sync else None
    ev_synced   = sync["events_synced"][()] if "events_synced" in sync else None


with h5py.File(out_path, "a") as f:
    # record which sync artifact was used
    f.attrs["sounds_sync_mode"] = mode

    if mode == "no_sync" or ev_synced is None:
        # nothing to merge; base remains canonical
        f.attrs["sounds_sync_mode"] = "no_sync"

    else:
        # --- read needed base datasets ---
        # base raw
        raw = f["raw"]
        sounds_logged = raw["sounds"][()]      # columns: [time, sound_id, ...?] (base assumes [t, id] at least)
        # processed timeline
        proc = f["processed"]
        timeline = proc["timeline"][()]        # columns include time at [:,0] and sound_ids at [:,6]
        times = timeline[:, 0]
        trials = proc["trial_idxs"][()]
        # parameters dict is stored as JSON in attrs
        parameters = json.loads(proc.attrs["parameters"])

        # --- ev_synced expected format ---
        # Your sync.py returns `ev_synced` with same columns as sounds.csv (it copies sounds_csv and shifts times),
        # so column 0 is sound time, column 1 is sound id.
        if ev_synced.ndim != 2 or ev_synced.shape[1] < 2:
            raise ValueError(f"events_synced has unexpected shape: {ev_synced.shape}")

        # store raw/sounds_ephys (like original pack.py idea)
        if "sounds_ephys" in raw:
            del raw["sounds_ephys"]
        ds = raw.create_dataset("sounds_ephys", data=ev_synced)
        ds.attrs["headers"] = raw["sounds"].attrs.get("headers", "")

        # --- rebuild processed/sound_events based on ev_synced ---
        sound_events = np.zeros((len(ev_synced), 3), dtype=float)
        for i in range(len(ev_synced)):
            t = float(ev_synced[i, 0])
            sid = float(ev_synced[i, 1])
            idx = nearest_tl_index(times, t)
            sound_events[i] = (t, sid, idx)

        # Remove events that map beyond timeline (defensive)
        sound_events = sound_events[sound_events[:, 2] < len(times)]

        if "sound_events" in proc:
            del proc["sound_events"]
        ds = proc.create_dataset("sound_events", data=sound_events)
        ds.attrs["headers"] = "sound_time, sound_id, timeline_idx"

        # --- update timeline sound_ids column only ---
        # base pack uses "sound_tl" where each timepoint holds the most recent played sound id
        # We'll recompute it from ev_synced.
        sound_tl = np.zeros(len(times), dtype=float)

        # sort by time (should already be sorted)
        ev_sorted = ev_synced[np.argsort(ev_synced[:, 0])]
        cur = 0
        cur_id = ev_sorted[0, 1] if len(ev_sorted) else 0
        for i, t in enumerate(times):
            while (cur + 1) < len(ev_sorted) and t >= ev_sorted[cur + 1, 0]:
                cur += 1
                cur_id = ev_sorted[cur, 1]
            sound_tl[i] = cur_id

        timeline[:, 6] = sound_tl  # 7th column in your base timeline headers

        # write back timeline
        del proc["timeline"]
        ds = proc.create_dataset("timeline", data=timeline)
        ds.attrs["headers"] = "time, x, y, speed, hd, trial_no, sound_ids, x_raw, y_raw"

        # --- rebuild target_matrix / distractor_matrix from updated sound_events ---
        tgt_matrix = build_tgt_matrix(sound_events, trials)
        if "target_matrix" in proc:
            del proc["target_matrix"]
        ds = proc.create_dataset("target_matrix", data=tgt_matrix)
        ds.attrs["headers"] = "sound_idx_start, sound_idx_end, tl_idx_start, tl_idx_end, result"

        if parameters.get("experiment", {}).get("distractor_fail", False):
            dis_matrix = build_dis_matrix(sound_events, trials, parameters)
            if "distractor_matrix" in proc:
                del proc["distractor_matrix"]
            ds = proc.create_dataset("distractor_matrix", data=dis_matrix)
            ds.attrs["headers"] = "sound_idx_start, sound_idx_end, tl_idx_start, tl_idx_end, result"

        # --- annotate correction ---
        proc.attrs["sound_time_correction"] = json.dumps(
            {
                "mode": mode,
                "source": os.path.basename(sync_path),
                "requested_sync": requested,
                "requested_type": requested_type,
                "used_ephys": use_ephys,
            }
        )
