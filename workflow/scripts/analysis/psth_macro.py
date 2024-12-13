
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
from utils.neurosuite import get_unit_names_sorted


with h5py.File(snakemake.input[0], 'r') as f:
    tl = np.array(f['processed']['timeline'])
    sound_events = np.array(f['processed']['sound_events'])
    tgt_matrix = np.array(f['processed']['target_matrix'])
    trials = np.array(f['processed']['trial_idxs'])
    cfg = json.loads(f['processed'].attrs['parameters'])

spike_times = {}
with h5py.File(snakemake.input[1], 'r') as f:
    unit_names = get_unit_names_sorted([name for name in f])
    for unit_name in f:
        spike_times[unit_name] = np.array(f[unit_name]['spike_times'])

tgt_dur = cfg['experiment']['target_duration']

# timeline idxs of success / miss target entrances
idxs_succ = np.where(tgt_matrix[:, 4] == 1)[0]
idxs_miss = np.where(tgt_matrix[:, 4] == 0)[0]

# noise offset is a silence without reward
noise_offset_idxs = []
for i in range(len(tl)):
    if tl[i, 6] == -1 and tl[i+1, 6] == 0:
        noise_offset_idxs.append(i+1)
noise_offset_times = tl[np.array(noise_offset_idxs)][:, 0]

# combinations to plot
dst_files = [x for x in snakemake.output]
macro_times = [
    [ tl[tgt_matrix[idxs_succ][:, 2]][:, 0], tl[tgt_matrix[idxs_miss][:, 2]][:, 0] ],  # target onset
    [ tl[tgt_matrix[idxs_succ][:, 3]][:, 0], tl[tgt_matrix[idxs_miss][:, 3]][:, 0] ],  # target offset
    [ tl[trials[:, 0].astype(np.int32)][:, 0] ],  # trial onset
    [ noise_offset_times ]  # noise offset
]
hw_bc = [[7, 51], [7, 51], [6, 49], [12, 49]]
spans = [(0, tgt_dur), (-tgt_dur, 0), (0, hw_bc[2][0]), (-10, 0)]
colors = [('green', 'black'), ('green', 'black'), ('tab:blue',), ('gray',)]
alphas = [(0.8, 0.5), (0.8, 0.5), (0.9,), (0.8,)]
labels = [('success', 'miss'), ('success', 'miss'), ('bgr',), ('noise',)]
titles = ['Target onset', 'Target offset', 'Trial onset', 'Noise offset']

for j, times in enumerate(macro_times):
    rows = int(np.ceil(len(unit_names)/3))
    fig = plt.figure(figsize=(15, rows*4))

    for i, unit_name in enumerate(unit_names):
        ax = fig.add_subplot(rows, 3, i+1)
        bins, counts_tgt_A = get_spike_counts(spike_times[unit_name], times[0], hw=hw_bc[j][0], bin_count=hw_bc[j][1])
        ax.hist(bins[:-1], bins=bins, weights=counts_tgt_A, edgecolor='black', color=colors[j][0], alpha=alphas[j][0], label=labels[j][0])
        if len(times) > 1:
            bins, counts_tgt_B = get_spike_counts(spike_times[unit_name], times[1], hw=hw_bc[j][0], bin_count=hw_bc[j][1])
            ax.hist(bins[:-1], bins=bins, weights=counts_tgt_B, edgecolor='black', color=colors[j][1], alpha=alphas[j][1], label=labels[j][1])
        ax.axvline(0, color='black', ls='--')
        ax.set_title("%s, %s" % (unit_name, titles[j]), fontsize=14)
        ax.axvspan(spans[j][0], spans[j][1], alpha=0.3, color='gray')
        ax.legend(loc='upper right', prop={'size': 10})
        if i % 3 == 0:
            ax.set_ylabel("Firing Rate, Hz", fontsize=14)

    fig.savefig(dst_files[j])