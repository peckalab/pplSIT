import os
import h5py
import json
import pywt
import numpy as np
from scipy import signal


area = snakemake.config['lfp']['aep_ITPC']['area']
metric = snakemake.config['lfp']['aep_ITPC']['metric']
fs = snakemake.config['lfp']['target_rate']

# AEPs
with h5py.File(snakemake.input[0], 'r') as f:
    lfp_trials = np.array(f[area][metric])  # (n_trials, n_samples)

# reading state indices
state_ids    = ['idxs_tgt_sta_succ', 'idxs_bgr_sta', 'idxs_bgr_run', 'idxs_sil_sta', 'idxs_sil_run']
state_labels = ['Target', 'Background (sta)', 'Background (run)', 'No stimulus (sta)', 'No stimulus (run)']
state_colors = ['tab:orange', 'tab:blue', 'red', 'grey', 'black']
with h5py.File(snakemake.input[1], 'r') as f:
    state_idxs = {}
    for idxs_name in state_ids:
        state_idxs[idxs_name] = np.array(f[idxs_name])


# === INPUTS ===
n_trials, n_samples = lfp_trials.shape
times = np.arange(lfp_trials.shape[1]) / fs * 1000  # in ms

# Define frequency range
frequencies = np.linspace(1, 79, 40)  # 1–80 Hz, 40 steps (2 Hz steps)
wavelet = 'cmor1.5-1.0'  # Complex Morlet: balance between time and frequency resolution

# Compute corresponding scales for desired frequencies
scales = pywt.central_frequency(wavelet) * fs / frequencies

collected = {}
for k, (state_name, idxs_state) in enumerate(state_idxs.items()):
    # === COMPUTE CWT ===
    # Returns: (n_trials, n_freqs, n_samples)
    cwt_all = np.array([
        pywt.cwt(trial, scales, wavelet, sampling_period=1/fs)[0]
        for trial in lfp_trials[idxs_state]
    ])
    
    # === ITPC ===
    # Phase normalization
    phases_real = cwt_all / np.abs(cwt_all)
    itpc_real = np.abs(np.mean(phases_real, axis=0))  # shape: (n_freqs, n_samples)
    
    # === POWER ===
    power_real = np.mean(np.abs(cwt_all) ** 2, axis=0)  # shape: (n_freqs, n_samples)
    
    # === Special: power for a DIP point at 50-70ms ===
    power_all = np.abs(cwt_all) ** 2  # shape: (n_trials, n_freqs, n_samples)

    freq_mask = (frequencies >= 1) & (frequencies <= 80)
    time_mask = (times >= 50) & (times <= 70)
   
    power_subset = power_all[:, freq_mask, :][:, :, time_mask]  # shape: (n_trials, selected_freqs, selected_times)
    power_avg_per_trial = power_subset.mean(axis=(0, 2))
    
    # === Compute trial-shuffled null ===
    n_shuffles = 100
    itpc_null = np.zeros((n_shuffles, len(frequencies), n_samples))
    
    for i in range(n_shuffles):
        shifted_phases = np.empty_like(phases_real)
    
        for trial in range(phases_real.shape[0]):
            shift = np.random.randint(0, n_samples)  # random circular shift
            shifted_phases[trial] = np.roll(phases_real[trial], shift=shift, axis=1)
    
        itpc_null[i] = np.abs(np.mean(shifted_phases, axis=0))
    
    # === STEP 3: Compute Z-score or p-value
    itpc_shuf_mean = np.mean(itpc_null, axis=0)
    itpc_shuf_std = np.std(itpc_null, axis=0)
    itpc_shuf_z = (itpc_real - itpc_shuf_mean) / (itpc_shuf_std + 1e-5)

    p_val = np.mean(itpc_null >= itpc_real[None, :, :], axis=0)
    significant_mask = p_val < 0.05

    itpc_threshold = np.percentile(itpc_null, 95, axis=0)  # shape: (freqs, times)
    non_phaselocked_mask = itpc_real < itpc_threshold

    # Optional: threshold significance at Z > 2 (roughly p < 0.05, two-sided)
    #significant_mask = itpc_z > 2

    # collect
    collected[state_name] = {
        'times': times,
        'frequencies': frequencies,
        'power_real': power_real,
        'itpc_real': itpc_real,
        'itpc_shuf_mean': itpc_shuf_mean,
        'itpc_shuf_std': itpc_shuf_std,
        'itpc_shuf_z': itpc_shuf_z,
        'significant_mask': significant_mask,
        'non_phaselocked_mask': non_phaselocked_mask,
        'dip_power': power_avg_per_trial
    }

with h5py.File(snakemake.output[0], 'w') as f:
    for state_name, idxs_state in state_idxs.items():
        grp = f.create_group(state_name)
        for k, ds in collected[state_name].items():
            grp.create_dataset(k, data=ds)
