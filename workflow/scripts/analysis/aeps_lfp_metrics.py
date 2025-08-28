import h5py, os, sys, json
import numpy as np



# reading data
with h5py.File(snakemake.input[0], 'r') as f:
    sound_events = np.array(f['processed']['sound_events'])
    tl = np.array(f['processed']['timeline'])
    cfg = json.loads(f['processed'].attrs['parameters'])

with h5py.File(snakemake.input[1], 'r') as f:
    lfp = np.array(f['lfp'])

# Set artifacts to 0. Ideally exclude these periods, but
# because it depends on the channel it's too complex
artifact_idxs_all = []
for i in range(len(lfp)):
    artf_idxs = np.where(np.abs(lfp[i]) > 4*lfp[i].std())[0]
    lfp[i][artf_idxs] = 0
    artifact_idxs_all.append(artf_idxs)

# time x channels
lfp = lfp.T

# baselines for each channel
with h5py.File(snakemake.input[2], 'r') as f:
    baselines = np.array(f['lfp_base'])
baseline_stds = baselines[:, 1]

fs = snakemake.config['lfp']['target_rate']  # Hz
pre_ms = snakemake.config['lfp']['aep_pre'] * 1000
post_ms = snakemake.config['lfp']['aep_dur'] * 1000
signal_window_ms = snakemake.config['lfp']['aep_sig'] * 1000
stim_times = sound_events[:, 0] * fs

n_ch = lfp.shape[1]
pre_samp = int(pre_ms * fs / 1000)
post_samp = int(post_ms * fs / 1000)
signal_win = int(signal_window_ms * fs / 1000)
win_len = pre_samp + post_samp

snr_trials = []
amp_trials = []
base_trials = []

for t in stim_times:
    t = int(t)
    if t - pre_samp >= 0 and t + post_samp < len(lfp):
        segment = lfp[t - pre_samp : t + post_samp, :]  # shape: time x channels

        baseline = segment[0:pre_samp, :]
        response = segment[pre_samp:pre_samp+signal_win, :]

        # Compute per-channel signal and noise
        base_amp = np.ptp(baseline, axis=0)
        signal = np.ptp(response, axis=0)     # peak-to-peak response
        #noise = np.std(baseline, axis=0)      # std of baseline pre-stim
        noise = baseline_stds
        
        snr = signal / (noise + 1e-6)         # avoid division by zero
        snr_trials.append(snr)
        amp_trials.append(signal)
        base_trials.append(base_amp)

snr_trials = np.stack(snr_trials, axis=0)  # shape: trials x channels
amp_trials = np.stack(amp_trials, axis=0)  # shape: trials x channels
base_trials = np.stack(base_trials, axis=0)  # shape: trials x channels

snrs  = np.mean(snr_trials, axis=0)    # mean SNR across trials → per channel
amps  = np.mean(amp_trials, axis=0)    # mean amplitude across trials → per channel
bases = np.mean(base_trials, axis=0)   # mean baseline across trials → per channel

with h5py.File(snakemake.output[0], 'w') as f:
    base_mx = np.vstack([snrs, amps, bases]).T
    f.create_dataset('aeps_lfp_metrics', data=base_mx)