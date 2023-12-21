import os, sys
import h5py
import json
import numpy as np
import matplotlib.pyplot as plt

# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.psth import get_spike_counts


with h5py.File(snakemake.input[0], 'r') as f:
    tl = np.array(f['processed']['timeline'])
    sound_events = np.array(f['processed']['sound_events'])
    cfg = json.loads(f['processed'].attrs['parameters'])

spike_times = {}
with h5py.File(snakemake.input[1], 'r') as f:
    unit_names = [name for name in f]
    for unit_name in f:
        spike_times[unit_name] = np.array(f[unit_name]['spike_times'])

bgr_dur    = cfg['sound']['sounds']['background']['duration']  # in seconds
dst_files  = [snakemake.output[0], snakemake.output[1]]
stim_combs = [(1, 2), (0, 1)]  # stimulus combinations
clr_combs  = [('tab:blue', 'tab:orange'), ('gray', 'tab:blue')]
titles     = [('bgr', 'tgt'), ('sil', 'bgr')]
hw = snakemake.config['psth']['micro']['latency']
bc = snakemake.config['psth']['micro']['bin_count']

for i, dst_file in enumerate(dst_files):
    psth_times_A = sound_events[sound_events[:, 1] == stim_combs[i][0]][:, 0]
    psth_times_B = sound_events[sound_events[:, 1] == stim_combs[i][1]][:, 0]

    cols = 3
    rows = int(np.ceil(len(unit_names)/3))
    fig = plt.figure(figsize=(15, rows*4))

    for j, unit_name in enumerate(unit_names):
        bins, counts_A = get_spike_counts(spike_times[unit_name], psth_times_A, hw, bc)
        bins, counts_B = get_spike_counts(spike_times[unit_name], psth_times_B, hw, bc)
        
        ax = fig.add_subplot(rows, cols, j+1)
        
        ax.hist(bins[:-1], bins=bins, weights=counts_A, edgecolor='black', color=clr_combs[i][0], alpha=0.7, label=titles[i][0])
        ax.hist(bins[:-1], bins=bins, weights=counts_B, edgecolor='black', color=clr_combs[i][1], alpha=0.7, label=titles[i][1])
        ax.axvline(0, color='black', ls='--')
        ax.axvspan(0, bgr_dur, alpha=0.3, color='gray')
        ax.axvspan(0 - hw, 0 - hw + bgr_dur, alpha=0.3, color='gray')
        ax.set_title(unit_name, fontsize=14)
        ax.legend(loc='lower right', prop={'size': 10})
        ax.set_xlim(-hw, hw)
        if j % 3 == 0:
            ax.set_ylabel("Firing Rate, Hz", fontsize=14)

    fig.savefig(dst_file)