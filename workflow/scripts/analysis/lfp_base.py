import h5py, os, sys, json
import numpy as np

# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.lfp import compute_bootstrapped_baseline


# reading data
with h5py.File(snakemake.input[0], 'r') as f:
    sound_events = np.array(f['processed']['sound_events'])
    tl = np.array(f['processed']['timeline'])
    cfg = json.loads(f['processed'].attrs['parameters'])
    trials = np.array(f['processed']['trial_idxs'])

with h5py.File(snakemake.input[1], 'r') as f:
    lfp = np.array(f['lfp'])

# clip artifacts
# for i in range(len(lfp)):
#     std = lfp[i].std()
#     lfp[i] = np.clip(lfp[i], a_min=-4*std, a_max=4*std)

# or set artifacts to 0. Ideally exclude these periods, but
# because it depends on the channel it's too complex
artifact_idxs_all = []
for i in range(len(lfp)):
    artf_idxs = np.where(np.abs(lfp[i]) > 4*lfp[i].std())[0]
    lfp[i][artf_idxs] = 0
    artifact_idxs_all.append(artf_idxs)

# time x channels
lfp = lfp.T

# Define thresholds
stationary_thresh = snakemake.config['lfp']['baseline']['stationary_thresh']  # m/s
min_duration_ms   = snakemake.config['lfp']['baseline']['min_duration_ms']  # only include segments longer than this
fs_lfp            = snakemake.config['lfp']['target_rate']  # Hz
fs_speed          = int(1/round(np.diff(tl[:, 0]).mean(), 6))  # Hz

# ITI-periods
ITI_start_idxs = trials[:, 1][:-1].astype(np.int32)
ITI_start_times = tl[ITI_start_idxs][:, 0]

ITI_end_idxs = trials[:, 0][1:].astype(np.int32)
ITI_end_times = tl[ITI_end_idxs][:, 0]

ITI_end_times = ITI_end_times[:len(ITI_start_times)]  # adjust in case last trial was not finished

inter_trial_period_times = np.vstack([ITI_start_times, ITI_end_times]).T
inter_trial_periods = (inter_trial_period_times * fs_lfp).astype(np.int32)  # in LFP samples

# speed
speed = tl[:, 3]

# Convert min duration to samples
min_duration_samples = int((min_duration_ms / 1000) * fs_lfp)

# Upsample speed to LFP resolution (linear interpolation)
from scipy.interpolate import interp1d

time_speed = np.linspace(0, len(speed)/fs_speed, len(speed))
time_lfp = np.linspace(0, len(lfp)/fs_lfp, len(lfp))

interp_func = interp1d(time_speed, speed, kind='linear', fill_value="extrapolate")
speed_upsampled = interp_func(time_lfp)  # same length as LFP time

# Find stationary and running segments during inter-trial periods
stationary_segments = []
running_segments = []

for start, end in inter_trial_periods:
    seg_speed = speed_upsampled[start:end]
    seg_lfp = lfp[start:end]

    # Find continuous blocks of stationary and running
    stationary_mask = seg_speed < stationary_thresh
    running_mask = seg_speed >= stationary_thresh

    def extract_segments(mask):
        segments = []
        in_segment = False
        seg_start = None
        for i, val in enumerate(mask):
            if val and not in_segment:
                in_segment = True
                seg_start = i
            elif not val and in_segment:
                in_segment = False
                if i - seg_start >= min_duration_samples:
                    segments.append((seg_start, i))
        if in_segment and (len(mask) - seg_start >= min_duration_samples):
            segments.append((seg_start, len(mask)))
        return segments

    stat_segs = extract_segments(stationary_mask)
    run_segs = extract_segments(running_mask)

    for s, e in stat_segs:
        stationary_segments.append(seg_lfp[s:e])
    for s, e in run_segs:
        running_segments.append(seg_lfp[s:e])

baseline_mean, baseline_std, baseline_ci_lower, baseline_ci_upper = compute_bootstrapped_baseline(stationary_segments, running_segments)

with h5py.File(snakemake.output[0], 'w') as f:
    base_mx = np.vstack([baseline_mean, baseline_std, baseline_ci_lower, baseline_ci_upper]).T
    f.create_dataset('lfp_base', data=base_mx)