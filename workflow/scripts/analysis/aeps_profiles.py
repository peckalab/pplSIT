import os, sys
import h5py
import numpy as np
import matplotlib.pyplot as plt

# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.aeps import outlier_lims, AEP_metrics_lims

# --------------------------

# load timeline / events
with h5py.File(snakemake.input[0], 'r') as f:
    tl           = np.array(f['processed']['timeline'])
    tgt_mx       = np.array(f['processed']['target_matrix'])
    sound_events = np.array(f['processed']['sound_events'])

# load AEPs    
aeps = {}
with h5py.File(snakemake.input[1], 'r') as f:
    for area in f:
        ds_name = [x for x in f[area] if x.find('filt') > 0][0]
        aeps[area] = np.array(f[area][ds_name])
        
# remove outliers - legacy for filtered LFP / AEPs
# for area, aeps_mx in aeps.items():
#     aeps_no_out = aeps_mx.copy()
#     aeps_mx[aeps_mx > outlier_lims[area]]  = outlier_lims[area]
#     aeps_mx[aeps_mx < -outlier_lims[area]] = -outlier_lims[area]
#     aeps[area] = aeps_no_out

# --------------------------

colors = [
    ['tab:blue', 'tab:orange'],
    ['tab:blue', 'green'],
    ['tab:blue', 'black'],
    ['tab:orange', 'black'],
    ['tab:blue', 'gray'],
]

labels = [
    ['BGR', 'TGT'],
    ['BGR', '1st success'],
    ['BGR', '1st reward'],
    ['TGT', '1st reward'],
    ['BGR', 'Silence'],
]

rows = len(labels)
cols = len(aeps.keys())
fig, axes = plt.subplots(rows, cols, figsize=(6*cols, 4.5*rows))

for k, (area, aeps_mx) in enumerate(aeps.items()):
    aeps_bgr_mean = aeps_mx[sound_events[:, 1] == 1].mean(axis=0)
    aeps_tgt_mean = aeps_mx[sound_events[:, 1] == 2].mean(axis=0)
    aeps_sil_mean = aeps_mx[sound_events[:, 1] == 0].mean(axis=0)
    aeps_1st_succ = aeps_mx[tgt_mx[tgt_mx[:, 4] == 1][:, 0]].mean(axis=0)
    aeps_1st_rewd = aeps_mx[tgt_mx[tgt_mx[:, 4] == 1][:, 1] + 1].mean(axis=0)

    combs = [
        [aeps_bgr_mean, aeps_tgt_mean],
        [aeps_bgr_mean, aeps_1st_succ],
        [aeps_bgr_mean, aeps_1st_rewd],
        [aeps_tgt_mean, aeps_1st_rewd],
        [aeps_bgr_mean, aeps_sil_mean],
    ]

    for i, comb in enumerate(combs):
        if cols == 1:
            ax = axes[i]
        else:
            ax = axes[i][k]
        if i == 0:
            ax.set_title(area, fontsize=14)
        ax.plot(comb[0] * 0.2, color=colors[i][0], label=labels[i][0])  # 0.2 OpenEphys scaling factor
        ax.plot(comb[1] * 0.2, color=colors[i][1], label=labels[i][1])
        ax.axhline(0, color='black')
        ax.axvline(0, color='black')
        ax.legend(loc='upper right', prop={'size': 14})
        ax.set_xlabel('Time from stim. onset, ms', fontsize=14)
        ax.grid()
        if k == 0:
            ax.set_ylabel('LFP, uV', fontsize=14)
        for j, (key, value) in enumerate(AEP_metrics_lims[area].items()):
            ax.axvline(value[0], color='black', ls='--', lw=1)
            ax.axvline(value[1], color='black', ls='--', lw=1)

fig.savefig(snakemake.output[0])
