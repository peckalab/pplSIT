import h5py, os, sys, json
import numpy as np
import matplotlib.pyplot as plt

# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.lfp import detect_artifacts_amp, clean_mask, global_rms_trace, highpass_lfp

# parameters
fs           = snakemake.config['lfp']['target_rate']  # Hz
win_ms       = snakemake.config['lfp']['artifacts']['win_ms']  # sliding window size for artifact detection
step_ms      = snakemake.config['lfp']['artifacts']['step_ms']  # sliding window step size for artifact detection
z_thresh     = snakemake.config['lfp']['artifacts']['z_thresh']  # z-score threshold for artifact detection
min_dur_ms   = snakemake.config['lfp']['artifacts']['min_dur_ms']  
merge_gap_ms = snakemake.config['lfp']['artifacts']['merge_gap_ms']

# load RAW lfp data
with h5py.File(snakemake.input[0], 'r') as f:
    lfp = np.array(f['lfp']).T  # time x channels
lfp = highpass_lfp(lfp, fs, 2.0)  # high-pass filter at 2 Hz

# compute global amplitude trace
global_amp = global_rms_trace(lfp)  # 1D

mask_amp, _, _ = detect_artifacts_amp(global_amp, fs, win_ms=win_ms, step_ms=step_ms, z_thresh=z_thresh)
artifact_mask  = clean_mask(mask_amp, fs, min_dur_ms=min_dur_ms, merge_gap_ms=merge_gap_ms)

# save artifacts to HDF5
with h5py.File(snakemake.output[0], 'w') as f:
    f.create_dataset('artifact_mask', data=artifact_mask.astype(np.uint8))  # save as uint8 to save space
    f.attrs['fs'] = fs
    f.attrs['win_ms'] = win_ms
    f.attrs['step_ms'] = step_ms
    f.attrs['z_thresh'] = z_thresh
    f.attrs['min_dur_ms'] = min_dur_ms
    f.attrs['merge_gap_ms'] = merge_gap_ms


# plot for debugging
fig, ax = plt.subplots(1, 1, figsize=(12, 5))

t = np.arange(len(lfp)) / fs

ax.plot(t, global_amp, label='LFP') # plot global amplitude trace
ax.fill_between(t, global_amp.min(), global_amp.max(), where=artifact_mask, alpha=0.2, color='red', label='Detected Artifacts')
ax.set_xlabel('Time (s)')
ax.set_ylabel('Global RMS Amplitude')
ax.set_title('Detected Artifacts (red shaded areas)')
ax.legend()
fig.tight_layout()
fig.savefig(snakemake.output[1])