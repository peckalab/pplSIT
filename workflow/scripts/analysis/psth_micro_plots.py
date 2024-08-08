import os, sys
import h5py, json
import numpy as np
import matplotlib.pyplot as plt

# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.neurosuite import get_unit_names_sorted


# reading some configs
with h5py.File(snakemake.input[0], 'r') as f:
    cfg = json.loads(f['processed'].attrs['parameters'])
with h5py.File(snakemake.input[1], 'r') as f:
    units_to_plot = get_unit_names_sorted([unit for unit in f['BGR']])  # should always exist

# configuration
event_types = [0, 1, 2, -1]  # SIL, BGR, TGT, NOI - order matters
colors = {0: 'gray', 1: 'tab:blue', 2: 'tab:orange', -1: 'red'}
ev_names = {0: 'SIL', 1: 'BGR', 2: 'TGT', -1: 'NOI'}
bgr_dur    = cfg['sound']['sounds']['background']['duration']  # in seconds
stim_combs = [(1, 2), (0, 1)]  # stimulus combinations to plot
cols = 3
rows = int(np.ceil(len(units_to_plot)/cols))
colors = {0: 'gray', 1: 'tab:blue', 2: 'tab:orange', -1: 'red'}
ev_names = {0: 'SIL', 1: 'BGR', 2: 'TGT', -1: 'NOI'}


# line plot figures
for fig_id, stim_comb in enumerate(stim_combs):
    idx_ev_1 = stim_comb[0]
    idx_ev_2 = stim_comb[1]

    # figure / file for each stimulus combination
    fig = plt.figure(figsize=(4*cols, 4*rows))

    for i, unit_name in enumerate(units_to_plot):
        with h5py.File(snakemake.input[1], 'r') as f:
            curr_stats_1 = np.array(f[ev_names[idx_ev_1]][unit_name]['profile_stats'])  # bins, profile mean, std, perc 5, perc 95
            curr_stats_2 = np.array(f[ev_names[idx_ev_2]][unit_name]['profile_stats'])  # bins, profile mean, std, perc 5, perc 95
        
        ax = fig.add_subplot(rows, cols, i+1)
        for j, c_stats in enumerate([curr_stats_1, curr_stats_2]):
            ax.plot(c_stats[0], c_stats[1], alpha=0.95, color=colors[stim_comb[j]], lw=2, label=ev_names[stim_comb[j]])
            ax.fill_between(c_stats[0], c_stats[3], c_stats[4], color=colors[stim_comb[j]], alpha=0.4)
        
        ax.set_ylim(bottom=0)
        ax.axvline(0, color='black', ls='--')
        ax.set_title("%s" % unit_name, fontsize=14)
        ax.axvspan(0, bgr_dur, alpha=0.3, color='gray')
        ax.legend()
        ax.grid()
        if i % cols == 0:
            ax.set_ylabel("Firing Rate, Hz", fontsize=14)
            
    fig.tight_layout()
    fig.savefig(snakemake.output[fig_id])


# bar plot figures
for fig_id, stim_comb in enumerate(stim_combs):
    idx_ev_1 = stim_comb[0]
    idx_ev_2 = stim_comb[1]

    # figure / file for each stimulus combination
    fig = plt.figure(figsize=(4*cols, 4*rows))

    for i, unit_name in enumerate(units_to_plot):
        with h5py.File(snakemake.input[1], 'r') as f:
            curr_stats_1 = np.array(f[ev_names[idx_ev_1]][unit_name]['profile_stats'])  # bins, profile mean, std, perc 5, perc 95
            curr_stats_2 = np.array(f[ev_names[idx_ev_2]][unit_name]['profile_stats'])  # bins, profile mean, std, perc 5, perc 95

        bins = list(curr_stats_1[0]) + [curr_stats_1[0][-1] + (curr_stats_1[0][1] - curr_stats_1[0][0])]
        hw = np.abs(curr_stats_1[0][0])
        
        ax = fig.add_subplot(rows, cols, i+1)
        ax.hist(bins[:-1], bins=bins, weights=curr_stats_1[1], edgecolor='black', color=colors[stim_comb[0]], alpha=0.7, label=ev_names[stim_comb[0]])
        ax.hist(bins[:-1], bins=bins, weights=curr_stats_2[1], edgecolor='black', color=colors[stim_comb[1]], alpha=0.7, label=ev_names[stim_comb[1]])
        ax.axvline(0, color='black', ls='--')
        ax.axvspan(0, bgr_dur, alpha=0.3, color='gray')
        ax.axvspan(0 - hw, 0 - hw + bgr_dur, alpha=0.3, color='gray')
        ax.set_title(unit_name, fontsize=14)
        ax.legend(loc='lower right', prop={'size': 10})
        ax.set_xlim(-hw, hw)
        if j % 3 == 0:
            ax.set_ylabel("Firing Rate, Hz", fontsize=14)
            
    fig.tight_layout()
    fig.savefig(snakemake.output[fig_id + 2])