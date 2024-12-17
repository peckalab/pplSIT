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
# Compute stationary_during_sound
stationary_thresh = snakemake.config['AEPs']['speed_thresh']
stationary_during_sound = np.array([np.all(tl[int(tl_idx):int(sound_events[idx+1,2]),3]<stationary_thresh) if idx<len(sound_events[:,2])-1 else np.all(tl[int(tl_idx):,3]<stationary_thresh) for idx, tl_idx in enumerate(sound_events[:,2])])

# load AEPs    
aeps = {}
with h5py.File(snakemake.input[1], 'r') as f:
    for area in f:
        ds_name = [x for x in f[area] if x.find('filt') > 0][0]
        aeps[area] = np.array(f[area][ds_name])

aep_dur = snakemake.config['lfp']['aep_dur']  # AEP duration in sec
        
# remove outliers - legacy for filtered LFP / AEPs
# for area, aeps_mx in aeps.items():
#     aeps_no_out = aeps_mx.copy()
#     aeps_mx[aeps_mx > outlier_lims[area]]  = outlier_lims[area]
#     aeps_mx[aeps_mx < -outlier_lims[area]] = -outlier_lims[area]
#     aeps[area] = aeps_no_out

# --------------------------

colors = [
    ['tab:blue', 'tab:orange'],
    ['navy', 'orangered'],
    ['navy', 'cyan'],
    ['tab:blue', 'green'],
    ['tab:blue', 'black'],
    ['tab:orange', 'black'],
    ['tab:blue', 'gray'],
]

labels = [
    ['BGR', 'TGT'],
    ['BGR (sta)', 'TGT (sta)'],
    ['BGR (sta)', 'BGR (run)'],
    ['BGR', '1st success'],
    ['BGR', '1st reward'],
    ['TGT', '1st reward'],
    ['BGR', 'Silence'],
]

rows = len(labels)
cols = len(aeps.keys())
fig, axes = plt.subplots(rows, cols, figsize=(6*cols, 4.5*rows))

for k, (area, aeps_mx) in enumerate(aeps.items()):
    # means
    aeps_bgr_mean = aeps_mx[sound_events[:, 1] == 1].mean(axis=0)
    aeps_tgt_mean = aeps_mx[sound_events[:, 1] == 2].mean(axis=0)
    aeps_bgr_sta_mean = aeps_mx[(sound_events[:, 1] == 1) & (stationary_during_sound)].mean(axis=0)
    aeps_bgr_run_mean = aeps_mx[(sound_events[:, 1] == 1) & (~stationary_during_sound)].mean(axis=0)
    aeps_tgt_sta_mean = aeps_mx[(sound_events[:, 1] == 2) & (stationary_during_sound)].mean(axis=0)
    aeps_sil_mean = aeps_mx[sound_events[:, 1] == 0].mean(axis=0)
    aeps_1st_succ = aeps_mx[tgt_mx[tgt_mx[:, 4] == 1][:, 0]].mean(axis=0)
    aeps_1st_rewd = aeps_mx[tgt_mx[tgt_mx[:, 4] == 1][:, 1] + 1].mean(axis=0)

    # counts
    aeps_bgr_count = np.sum(sound_events[:, 1] == 1)
    aeps_tgt_count = np.sum(sound_events[:, 1] == 2)
    aeps_bgr_sta_count = np.sum((sound_events[:, 1] == 1) & (stationary_during_sound))
    aeps_bgr_run_count = np.sum((sound_events[:, 1] == 1) & (~stationary_during_sound))
    aeps_tgt_sta_count = np.sum((sound_events[:, 1] == 2) & (stationary_during_sound))
    aeps_sil_count = np.sum(sound_events[:, 1] == 0)
    aeps_1st_succ_count = len(tgt_mx[tgt_mx[:, 4] == 1][:, 0])
    aeps_1st_rewd_count = len(tgt_mx[tgt_mx[:, 4] == 1][:, 1] + 1)

    # STEs
    aeps_bgr_ste = aeps_mx[sound_events[:, 1] == 1].std(axis=0) / np.sqrt(aeps_bgr_count)
    aeps_tgt_ste = aeps_mx[sound_events[:, 1] == 2].std(axis=0) / np.sqrt(aeps_tgt_count)
    aeps_bgr_sta_ste = aeps_mx[(sound_events[:, 1] == 1) & (stationary_during_sound)].std(axis=0) / np.sqrt(aeps_bgr_sta_count)
    aeps_bgr_run_ste = aeps_mx[(sound_events[:, 1] == 1) & (~stationary_during_sound)].std(axis=0) / np.sqrt(aeps_bgr_run_count)
    aeps_tgt_sta_ste = aeps_mx[(sound_events[:, 1] == 2) & (stationary_during_sound)].std(axis=0) / np.sqrt(aeps_tgt_sta_count)
    aeps_sil_ste = aeps_mx[sound_events[:, 1] == 0].std(axis=0) / np.sqrt(aeps_sil_count)
    aeps_1st_succ_ste = aeps_mx[tgt_mx[tgt_mx[:, 4] == 1][:, 0]].std(axis=0) / np.sqrt(aeps_1st_succ_count)
    aeps_1st_rewd_ste = aeps_mx[tgt_mx[tgt_mx[:, 4] == 1][:, 1] + 1].std(axis=0) / np.sqrt(aeps_1st_rewd_count)


    combs = [
        [aeps_bgr_mean, aeps_tgt_mean],
        [aeps_bgr_sta_mean, aeps_tgt_sta_mean],
        [aeps_bgr_sta_mean, aeps_bgr_run_mean],
        [aeps_bgr_mean, aeps_1st_succ],
        [aeps_bgr_mean, aeps_1st_rewd],
        [aeps_tgt_mean, aeps_1st_rewd],
        [aeps_bgr_mean, aeps_sil_mean],
    ]

    combs_counts = [
        [aeps_bgr_count, aeps_tgt_count],
        [aeps_bgr_sta_count, aeps_tgt_sta_count],
        [aeps_bgr_sta_count, aeps_bgr_run_count],
        [aeps_bgr_count, aeps_1st_succ_count],
        [aeps_bgr_count, aeps_1st_rewd_count],
        [aeps_tgt_count, aeps_1st_rewd_count],
        [aeps_bgr_count, aeps_sil_count],
    ]

    combs_ste = [
        [aeps_bgr_ste, aeps_tgt_ste],
        [aeps_bgr_sta_ste, aeps_tgt_sta_ste],
        [aeps_bgr_sta_ste, aeps_bgr_run_ste],
        [aeps_bgr_ste, aeps_1st_succ_ste],
        [aeps_bgr_ste, aeps_1st_rewd_ste],
        [aeps_tgt_ste, aeps_1st_rewd_ste],
        [aeps_bgr_ste, aeps_sil_ste],
    ]

    for i, comb in enumerate(combs):
        if cols == 1:
            ax = axes[i]
        else:
            ax = axes[i][k]
        if i == 0:
            ax.set_title(area, fontsize=14)
        

        # shade with STE
        ax.fill_between(np.linspace(0,aep_dur*1000,len(combs_ste[i][0])),0.2*(comb[0] - combs_ste[i][0]), 0.2*(comb[0]+combs_ste[i][0]), color=colors[i][0], alpha=0.3)
        # line plot of the mean
        ax.plot(comb[0] * 0.2, color=colors[i][0], label=f'{labels[i][0]} ({combs_counts[i][0]})')  # 0.2 OpenEphys scaling factor
        
        # shade with STE
        ax.fill_between(np.linspace(0,aep_dur*1000,len(combs_ste[i][1])),0.2*(comb[1] - combs_ste[i][1]), 0.2*(comb[1]+combs_ste[i][1]), color=colors[i][1], alpha=0.3)
        # line plot of the mean
        ax.plot(comb[1] * 0.2, color=colors[i][1], label=f'{labels[i][1]} ({combs_counts[i][1]})')
        
        ax.axhline(0, color='black')
        ax.axvline(0, color='black')
        ax.legend(loc='upper right', prop={'size': 14})
        ax.set_xlabel('Time from stim. onset, ms', fontsize=14)
        ax.grid()
        if k == 0:
            ax.set_ylabel(r'LFP, $\mu$V', fontsize=14)
        for j, (key, value) in enumerate(AEP_metrics_lims[area].items()):
            ax.axvline(value[0], color='black', ls='--', lw=1)
            ax.axvline(value[1], color='black', ls='--', lw=1)

fig.savefig(snakemake.output[0])
