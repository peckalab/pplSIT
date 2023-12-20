import h5py
import time
import os, json
import numpy as np
from scipy import signal


def head_direction(tl, hd_update_speed=0.04):
    width = 200  # 100 points ~= 1 sec with at 100Hz
    kernel = signal.gaussian(width, std=(width) / 7.2)

    x_smooth = np.convolve(tl[:, 1], kernel, 'same') / kernel.sum()
    y_smooth = np.convolve(tl[:, 2], kernel, 'same') / kernel.sum()

    diff_x = np.diff(x_smooth, axis=0)
    diff_y = np.diff(y_smooth, axis=0)

    hd = -np.arctan2(diff_y, diff_x)
    hd = np.concatenate([np.array([hd[0]]), hd])  # make the same length as timeline
    
    # reset idle periods
    idle_periods = []
    idle_idxs = np.where(tl[:, 3] < hd_update_speed)[0]

    crit = np.where(np.diff(idle_idxs) > 1)[0]

    idle_periods.append( (idle_idxs[0], idle_idxs[crit[0]]) )  # first idle period
    for i, point in enumerate(crit[:-1]):
        idx_start = idle_idxs[crit[i] + 1]
        idx_end = idle_idxs[crit[i+1]]
        idle_periods.append( (idx_start, idx_end) )

    idle_periods = np.array(idle_periods)
    
    for (i1, i2) in idle_periods:
        hd[i1:i2] = hd[i1-1]
        
    return hd


def build_tgt_matrix(sound_events, trials):
    # compute timeline / sound indices of entrances / exits to the target
    tgt_start_idxs = []
    tgt_end_idxs = []
    
    for i, event in enumerate(sound_events[:-1]):
        if sound_events[i][1] < 2 and sound_events[i+1][1] == 2:
            tgt_start_idxs.append(i+1)
        if sound_events[i][1] == 2 and sound_events[i+1][1] < 2:
            tgt_end_idxs.append(i)
    
    # ignore first/last target if not ended
    if tgt_start_idxs[-1] > tgt_end_idxs[-1]:
        tgt_start_idxs = tgt_start_idxs[:-1]
    if tgt_end_idxs[0] < tgt_start_idxs[0]:
        tgt_end_idxs = tgt_end_idxs[1:]
    tgt_start_idxs = np.array(tgt_start_idxs)
    tgt_end_idxs   = np.array(tgt_end_idxs)
    
    # successful / missed
    tgt_results = np.zeros(len(tgt_start_idxs))
    for idx_tl_success_end in trials[trials[:, 5] == 1][:, 1]:
        idx_succ = np.abs(sound_events[tgt_end_idxs][:, 2] - idx_tl_success_end).argmin()
        tgt_results[idx_succ] = 1
        
    # tl_idx_start, tl_idx_end, aep_idx_start, aer_idx_end, success / miss
    return np.column_stack([
        tgt_start_idxs,
        tgt_end_idxs,
        sound_events[tgt_start_idxs][:, 2],
        sound_events[tgt_end_idxs][:, 2],
        tgt_results
    ]).astype(np.int32)


def pack(pos_file, ev_file, snd_file, cfg_file, man_file, dst_file, drift_coeff=0.000025):  
    """
    Pack independent raw session datasets into a single HDF5 file.

    args:
        pos_file    - path to positions file
        ev_file     - path to events file
        snd_file    - path to sounds file
        cfg_file    - path to config file
        dst_file    - destination HDF5 file
        drift_coeff - 15 ms per 10 minutes default sound drift
    
    File has the following structure:
    
    /raw
        /positions      - raw positions from .csv
        /events         - raw events from .csv
        /sounds         - raw sounds from .csv
        /islands        - raw island infos from .csv (if exists)
    /processed
        /timeline       - matrix of [time, x, y, speed, HD, trial_no, sound_id] sampled at 100Hz,
                          data is smoothed using gaussian kernels,
                          inter-trial intervals have trial_no = 0
        /trial_idxs     - matrix of trial indices to timeline
        /sound_events     - matrix of sound times and indices to timeline
        
    each dataset has an attribute 'headers' with the description of columns.
    """
    with open(os.path.join(cfg_file)) as json_file:
        parameters = json.load(json_file)
    
    with h5py.File(dst_file, 'w') as f:  # overwrite mode

        # -------- save raw data ------------
        raw = f.create_group('raw')
        raw.attrs['parameters'] = json.dumps(parameters)

        ds_names = ['positions', 'events', 'sounds', 'islands']
        for i, f_path  in enumerate([pos_file, ev_file, snd_file]):
            ds_name = ds_names[i]
            if not os.path.exists(f_path):
                continue
                
            with open(f_path) as ff:
                headers = ff.readline()
            data = np.loadtxt(f_path, delimiter=',', skiprows=1)

            ds = raw.create_dataset(ds_name, data=data)
            ds.attrs['headers'] = headers
        
        # read raw data and normalize to session start
        events = np.array(f['raw']['events'])
        s_start, s_end = events[:, 0][0], events[:, 0][-1]
        events[:, 0] = events[:, 0] - s_start
        positions = np.array(f['raw']['positions'])
        positions[:, 0] = positions[:, 0] - s_start
        sounds = np.array(f['raw']['sounds'])
        sounds[:, 0] = sounds[:, 0] - s_start

        # squeeze - if session was interrupted, adjust times
        # to have a continuous timeline
        end_idxs = np.where(events[:, 5] == -1)[0]
        if len(end_idxs) > 1:
            # diffs in time beetween pauses
            deltas = [events[idx + 1][0] - events[idx][0] for idx in end_idxs[:-1]]

            for df, delta in zip(end_idxs, deltas):  # squeezing events
                events[df+1:][:, 0] = events[df+1:][:, 0] - delta

            end_idxs = np.where(np.diff(sounds[:, 0]) > 20)[0]  # squeezing sounds
            for df, delta in zip(end_idxs, deltas):
                sounds[df+1:][:, 0] = sounds[df+1:][:, 0] - delta

            end_idxs = np.where(np.diff(positions[:, 0]) > 20)[0]  # squeezing positions - more than 20? secs pauses
            for df, delta in zip(end_idxs, deltas):
                positions[df+1:][:, 0] = positions[df+1:][:, 0] - delta
            parameters['experiment']['timepoints'] = [positions[df+1][0] for df in end_idxs]  # update session parameters
            parameters['experiment']['session_duration'] = positions[-1][0]

        # -------- save processed ------------
        proc = f.create_group('processed')
        proc.attrs['parameters'] = json.dumps(parameters)

        # convert timeline to 100 Hz
        time_freq = 100  # at 100Hz
        s_start, s_end = events[0][0], events[-1][0]
        times = np.linspace(s_start, s_end, int((s_end - s_start) * time_freq))
        pos_at_freq = np.zeros((len(times), 3))

        curr_idx = 0
        for i, t in enumerate(times):
            if curr_idx < len(positions) - 1 and \
                np.abs(t - positions[:, 0][curr_idx]) > np.abs(t - positions[:, 0][curr_idx + 1]):
                curr_idx += 1
            pos_at_freq[i] = (t, positions[curr_idx][1], positions[curr_idx][2])

        # save trials
        t_count = len(np.unique(events[events[:, -1] != 0][:, -2]))
        trials = np.zeros((t_count, 6))
        for i in range(t_count):
            t_start_idx = (np.abs(pos_at_freq[:, 0] - events[2*i][0])).argmin()
            t_end_idx = (np.abs(pos_at_freq[:, 0] - events[2*i + 1][0])).argmin()
            state = 1 if events[2*i + 1][-1] == 1 else 0

            trials[i] = (t_start_idx, t_end_idx, events[2*i][1], events[2*i][2], events[2*i][3], state)

        trial_idxs = proc.create_dataset('trial_idxs', data=trials)
        trial_idxs.attrs['headers'] = 't_start_idx, t_end_idx, target_x, target_y, target_r, fail_or_success'

        # adjust for a drift and offset
        drift = s_end * drift_coeff
        with open(man_file, 'r') as json_file:
            offset = json.load(json_file)['ephys']['offset']
        sounds[:, 0] = sounds[:, 0] + np.arange(len(sounds)) * drift/len(sounds) + offset/1000.

        # save sounds
        sound_events = np.zeros((len(sounds), 3))
        left_idx = 0
        delta = 10**5
        for i in range(len(sounds)):
            while left_idx < len(pos_at_freq) and \
                    np.abs(sounds[i][0] - pos_at_freq[:, 0][left_idx]) < delta:
                delta = np.abs(sounds[i][0] - pos_at_freq[:, 0][left_idx])
                left_idx += 1

            sound_events[i] = (sounds[i][0], sounds[i][1], left_idx)
            delta = 10**5

        sound_events_ds = proc.create_dataset('sound_events', data=sound_events)
        sound_events_ds.attrs['headers'] = 'sound_time, sound_id, timeline_idx'

        # building timeline
        width = 50  # 100 points ~= 1 sec with at 100Hz
        kernel = signal.gaussian(width, std=(width) / 7.2)

        x_smooth = np.convolve(pos_at_freq[:, 1], kernel, 'same') / kernel.sum()
        y_smooth = np.convolve(pos_at_freq[:, 2], kernel, 'same') / kernel.sum()

        # speed
        dx = np.sqrt(np.square(np.diff(x_smooth)) + np.square(np.diff(y_smooth)))
        dt = np.diff(pos_at_freq[:, 0])
        speed = np.concatenate([dx/dt, [dx[-1]/dt[-1]]])

        # head direction - zeros for the moment
        temp_tl = np.column_stack([pos_at_freq[:, 0], x_smooth, y_smooth, speed])
        hd = head_direction(temp_tl)

        # trial numbers
        trials_data = np.zeros(len(temp_tl))

        for i, trial in enumerate(trials):
            idx1, idx2 = trial[0], trial[1]
            trials_data[int(idx1):int(idx2)] = i + 1

        # sounds played
        sound_tl = np.zeros(len(temp_tl))
        curr_sound_idx = 0
        for i in range(len(temp_tl)):
            if curr_sound_idx + 1 >= len(sounds):
                break

            if temp_tl[i][0] > sounds[curr_sound_idx][0]:
                curr_sound_idx += 1
            sound_tl[i] = sounds[curr_sound_idx][1]

        timeline = proc.create_dataset('timeline', data=np.column_stack(\
                     [pos_at_freq[:, 0], x_smooth, y_smooth, speed, hd, trials_data, sound_tl, x_smooth, y_smooth],
                   ))
        timeline.attrs['headers'] = 'time, x, y, speed, hd, trial_no, sound_ids, x_raw, y_raw'

        # create target matrix
        tgt_matrix = build_tgt_matrix(sound_events, trials)
        tgt_matrix_ds = proc.create_dataset('target_matrix', data=tgt_matrix)
        tgt_matrix_ds.attrs['headers'] = 'sound_idx_start, sound_idx_end, tl_idx_start, tl_idx_end, result'


# actual execution
pack(
    snakemake.input[0], 
    snakemake.input[1], 
    snakemake.input[2], 
    snakemake.input[3], 
    snakemake.input[4],
    snakemake.output[0]
)
