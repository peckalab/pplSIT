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
colors = {'SIL': 'gray', 'BGR': 'tab:blue', 'TGT': 'tab:orange', 'NOI': 'red', 'DI1': 'tab:green', 'idxs_bgr_sta': 'tab:blue', 
          'idxs_tgt_sta': 'tab:orange', 'idxs_di1_sta': 'tab:green', 'idxs_tgt_sta_succ': 'tab:orange', 
          'idxs_dis_fail': 'tab:green'}
ev_names = {'SIL': 'SIL', 'BGR': 'BGR', 'TGT': 'TGT', 'NOI': 'NOI', 'DI1': 'DI1', 'idxs_bgr_sta': 'BGR', 
          'idxs_tgt_sta': 'TGT', 'idxs_di1_sta': 'DI1', 'idxs_tgt_sta_succ': 'TGT', 
          'idxs_dis_fail': 'DIS'}
bgr_dur    = cfg['sound']['sounds']['background']['duration']  # in seconds

stim_combs = [('BGR', 'TGT'), ('SIL', 'BGR')]  # stimulus combinations to plot - these always exist
fig_directory = os.path.dirname(snakemake.output[0])
fig_names = ["psth_bgr_tgt_line.pdf", "psth_bgr_sil_line.pdf"]

# get available profiles
with h5py.File(snakemake.input[1], 'r') as f:
    available_conds = [x for x in f]  # should be like ['BGR', 'TGT', 'SIL', 'NOI']

if ('idxs_bgr_sta' in available_conds)&('idxs_tgt_sta' in available_conds):
    stim_combs.append(('idxs_bgr_sta', 'idxs_tgt_sta'))
    fig_names.append("psth_bgr_tgt_sta_line.pdf")
    if 'idxs_di1_sta' in available_conds:
        stim_combs.append(('idxs_bgr_sta', 'idxs_tgt_sta', 'idxs_di1_sta'))
        fig_names.append("psth_bgr_tgt_di1_sta_line.pdf")
if 'DI1' in available_conds:
    stim_combs.append(('BGR', 'TGT', 'DI1'))
    fig_names.append("psth_bgr_tgt_di1_line.pdf")
if ('idxs_tgt_sta_succ' in available_conds)&('idxs_dis_fail' in available_conds):
    stim_combs.append(('idxs_tgt_sta_succ', 'idxs_dis_fail'))
    fig_names.append("psth_tgt_succ_dis_fail_line.pdf")

cols = 3
rows = int(np.ceil(len(units_to_plot)/cols))
latency = cfg['sound']['latency']


# line plot figures
for fig_id, stim_comb in enumerate(stim_combs):

    # figure / file for each stimulus combination
    fig = plt.figure(figsize=(4*cols, 4*rows))
    
    for i, unit_name in enumerate(units_to_plot):
        curr_stats = []
        with h5py.File(snakemake.input[1], 'r') as f:
            for stim_i in range(len(stim_comb)):
                curr_stats.append(np.array(f[stim_comb[stim_i]][unit_name]['profile_stats']))  # bins, profile mean, std, perc 5, perc 95
        
        ax = fig.add_subplot(rows, cols, i+1)
        for j, c_stats in enumerate(curr_stats):
            bin_size = c_stats[0][1] - c_stats[0][0]
            ax.plot(c_stats[0] + bin_size, c_stats[1], alpha=0.95, color=colors[stim_comb[j]], lw=2, label=ev_names[stim_comb[j]])
            ax.fill_between(c_stats[0] + bin_size, c_stats[3], c_stats[4], color=colors[stim_comb[j]], alpha=0.4)
        ax.set_xlim(-latency, latency)
        ax.set_ylim(bottom=0)
        ax.axvline(0, color='black', ls='--')
        ax.set_title("%s" % unit_name, fontsize=14)
        ax.axvspan(0, bgr_dur, alpha=0.3, color='gray')
        ax.axvspan(-latency, -latency + bgr_dur, alpha=0.3, color='gray')
        ax.legend()
        ax.grid()
        if i % cols == 0:
            ax.set_ylabel("Firing Rate, Hz", fontsize=14)
            
    fig.tight_layout()
    fig.savefig(os.path.join(fig_directory, fig_names[fig_id]))
