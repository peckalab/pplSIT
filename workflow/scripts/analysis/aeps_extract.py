import os
import h5py
import json
import numpy as np
from scipy import signal



# selected channels for AEPs
with open(snakemake.input[0]) as json_file:
    channels_all = json.load(json_file)['AEPs']

with h5py.File(snakemake.input[1], 'r') as f:
    lfp = np.array(f['lfp'])  # channels x time

with h5py.File(snakemake.input[2], 'r') as f:
    sound_events = np.array(f['processed']['sound_events'])

# LFP baselines
with h5py.File(snakemake.input[3], 'r') as f:
    baselines = np.array(f['lfp_base'])
    baseline_per_channel = baselines[:, 0]
    baseline_stds = baselines[:, 1]

# channels response metrics - amplitude / SNRs
with h5py.File(snakemake.input[4], 'r') as f:
    metrics = np.array(f['aeps_lfp_metrics'])  # snrs, amplitudes, pre-stim baselines - per channel
    amplitude_weights = metrics[:, 1]
    
    # Normalize weights
    weights = amplitude_weights / (np.sum(amplitude_weights) + 1e-6)

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
nyquist = 0.5 * fs


# clean LFP
for i in range(len(lfp)):
    lfp_ch = lfp[i]  # raw LFP

    # first cut out outliers (is this correct?)
    lfp_filt = np.clip(lfp_ch, a_min=-4*lfp_ch.std(), a_max=4*lfp_ch.std())

    # try to filter out 50Hz
    low  = 47 / nyquist
    high = 53 / nyquist
    sos = signal.butter(10, [low, high], btype='bandstop', output='sos')
    lfp_filt = signal.sosfiltfilt(sos, lfp_filt)

    # try to filter < 1Hz slow staff
    h_cut = 1 / nyquist
    sos = signal.butter(10, h_cut, btype='highpass', output='sos')
    lfp_filt = signal.sosfiltfilt(sos, lfp_filt)

    lfp[i] = lfp_filt

# to time x channels
lfp = lfp.T


for area, channels in channels_all.items():  # selected_channels are also indices to channel metrics
    if not type(channels) == list:
        selected_channels = [channels]
    else:
        selected_channels = channels

    avg_across_channels = []
    weighted_avg = []
    zscored_avg = []
    raw_AEPs = []  # pulse, samples, channels

    # extract responses
    for t in stim_times:
        t = int(t)
        #if t - pre_samp >= 0 and t + post_samp < len(lfp):
        if t + post_samp < len(lfp):  # last pulse is the same as pre-last, if overflows
            seg = lfp[t : t + post_samp, selected_channels]  # shape: time x channels

        # Baseline subtraction
        seg_baselined = seg - baseline_per_channel[selected_channels][None, :]

        # (1) Raw average
        avg_across_channels.append(np.mean(seg_baselined, axis=1))  # time

        # (2) Weighted average
        weighted_avg.append(np.sum(seg_baselined * weights[selected_channels][None, :], axis=1))  # time

        # (3) Z-scored average (per channel)
        seg_z = seg_baselined / baseline_stds[None, selected_channels]
        zscored_avg.append(np.mean(seg_z, axis=1))  # time

        # (4) Raw AEPs from selected channels
        raw_AEPs.append(seg_baselined)

        # (4) Peak-to-peak amplitude in signal window per channel
        #ptp = np.ptp(seg_baselined[:signal_win, :], axis=0)
        #EV_matrix.append(ptp)

        # (5) Peak-to-peak Mean amplitude in sustained window per channel
        #late = np.mean(seg_baselined[signal_win:, :], axis=0)
        #SU_matrix.append(late)

    avg_across_channels = np.stack(avg_across_channels, axis=0)
    weighted_avg = np.stack(weighted_avg, axis=0)
    zscored_avg = np.stack(zscored_avg, axis=0)
    raw_AEPs = np.stack(raw_AEPs, axis=0)

    # save AEPs to file
    with h5py.File(snakemake.output[0], 'a') as f:
        if not area in f:
            f.create_group(area)
        if len(f[area]) > 0:
            for name in f[area]:
                del f[area][name]

        f[area].create_dataset('avg_across_channels', data=avg_across_channels)
        f[area].create_dataset('weighted_avg', data=weighted_avg)
        f[area].create_dataset('zscored_avg', data=zscored_avg)
        f[area].create_dataset('raw_AEPs', data=raw_AEPs)  # pulse, samples, channels