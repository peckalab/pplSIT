import os, sys
import h5py, json
import numpy as np
import matplotlib.pyplot as plt

# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.neurosuite import get_unit_names_sorted
from utils.psth import get_spike_counts
from utils.maths import pval2text


s_path  = os.path.dirname(snakemake.input[0])
session = os.path.basename(s_path)
sound_phase_lock_file = os.path.join(s_path, 'analysis', 'sound_phase_lock.h5')

# reading some configs
with h5py.File(snakemake.input[0], 'r') as f:
    cfg = json.loads(f['processed'].attrs['parameters'])
    sound_events = np.array(f['processed']['sound_events'])
    tl = np.array(f['processed']['timeline'])

spike_times = {}
depth = {}
with h5py.File(snakemake.input[1], 'r') as f:
    units_to_plot = get_unit_names_sorted([name for name in f])
    for unit_name in units_to_plot:
        spike_times[unit_name] = np.array(f[unit_name]['spike_times'])
        if 'kilosort_info' in f[unit_name]:
            # depth is in um
            depth[unit_name] = f[unit_name]['kilosort_info'][4]
        else:
            # depth is in um
            depth[unit_name] = np.array(f[unit_name]['anatomical_position'])[1]
    # ordered by shank, then depth
    units_to_plot = sorted(units_to_plot, key=lambda x: (int(x.split('-')[0]), depth[x]))
    # units_to_plot = sorted(spike_times.keys(), key=lambda x: depth[x])

# configuration
event_types = [0, 1, 2, -1]  # SIL, BGR, TGT, NOI - order matters
colors = {0: 'gray', 1: 'tab:blue', 2: 'tab:orange', -1: 'red'}
ev_names = {0: 'SIL', 1: 'BGR', 2: 'TGT', -1: 'NOI'}

if 'background' in cfg['sound']['sounds']:
    bgr_dur    = cfg['sound']['sounds']['background']['duration']  # in seconds
    passive = False
else: # if background is not present, it is a passive experiment
    passive = True
    idxs_ev_passive = {}
    bgr_dur    = cfg['sound']['sounds']['F1']['duration']
    sound_ev_types = [ev for ev in cfg['sound']['sounds'].keys() if ev != 'noise']
    for i, ev in enumerate(sound_ev_types):
        idxs_ev_passive[ev] = np.where(sound_events[:, 1] == i+1)[0]


#stim_combs = [(1, 2), (0, 1)]  # stimulus combinations to plot
cols = 3
rows = int(np.ceil(len(units_to_plot)/cols))
colors = {0: 'gray', 1: 'tab:blue', 2: 'tab:orange', -1: 'red'}
ev_names = {0: 'SIL', 1: 'BGR', 2: 'TGT', -1: 'NOI'}
#latency = cfg['sound']['latency']

hw = snakemake.config['psth']['micro']['latency']
bc = snakemake.config['psth']['micro']['bin_count']

speed_max = 0.04
speed_ev = tl[sound_events[:, 2].astype(np.int32)][:, 3]
idxs_sta_ev = np.where(speed_ev < speed_max)[0]
idxs_run_ev = np.where(speed_ev > speed_max)[0]
idxs_bgr_ev = np.where(sound_events[:, 1] == 1)[0]
idxs_sil_ev = np.where(sound_events[:, 1] == 0)[0]
idxs_tgt_ev = np.where(sound_events[:, 1] == 2)[0]
idxs_di1_ev = np.where(sound_events[:, 1] == 3)[0]
idxs_di2_ev = np.where(sound_events[:, 1] == 4)[0]  # for the moment limit to 2 distractors

stim_comb_idxs = [
    [idxs_bgr_ev, idxs_tgt_ev],  # standard BGR / TGT
    [idxs_sil_ev, idxs_bgr_ev],  # standard SIL / BGR
    [np.intersect1d(idxs_bgr_ev, idxs_sta_ev), np.intersect1d(idxs_tgt_ev, idxs_sta_ev)],  # BGR stationary / TGT stationary
    [np.intersect1d(idxs_bgr_ev, idxs_sta_ev), np.intersect1d(idxs_bgr_ev, idxs_run_ev)],  # BGR stationary / BGR run
    [np.intersect1d(idxs_sil_ev, idxs_sta_ev), np.intersect1d(idxs_sil_ev, idxs_run_ev)],  # SIL stationary / SIL run
]

label_combs = [
    ['bgr', 'tgt'],
    ['sil', 'bgr'],
    ['bgr_sta', 'tgt_sta'],
    ['bgr_sta', 'bgr_run'],
    ['sil_sta', 'sil_run'],
]

color_combs = [
    ['tab:blue', 'tab:orange'],
    ['grey', 'tab:blue'],
    ['tab:blue', 'tab:orange'],
    ['navy', 'tab:blue'],
    ['grey', 'tab:red'],
]

# TGT, BGR, SIL bar plot figures
for fig_id, stim_comb in enumerate(stim_comb_idxs):
    #idx_ev_1 = stim_comb[0]
    #idx_ev_2 = stim_comb[1]

    #idxs_ev_1 = np.where(sound_events[:, 1] == event_types[idx_ev_1])[0]
    #idxs_ev_2 = np.where(sound_events[:, 1] == event_types[idx_ev_2])[0]

    idxs_ev_1 = stim_comb[0]
    idxs_ev_2 = stim_comb[1]

    # figure / file for each stimulus combination
    fig = plt.figure(figsize=(4*cols, 4*rows))

    for i, unit_name in enumerate(units_to_plot):
        bins, psth1 = get_spike_counts(spike_times[unit_name], sound_events[idxs_ev_1][:, 0], hw=hw, bin_count=bc)
        bins, psth2 = get_spike_counts(spike_times[unit_name], sound_events[idxs_ev_2][:, 0], hw=hw, bin_count=bc)
        
        # if sound phase locking exists - plot with the label
        label1 = label_combs[fig_id][0]
        label2 = label_combs[fig_id][1]
        if sound_phase_lock_file:
            with h5py.File(sound_phase_lock_file, 'r') as snd_f:
                if label_combs[fig_id][0] in snd_f:
                    MRL = np.array(snd_f[label_combs[fig_id][0]][unit_name]['MRL_real'])
                    pv  = np.array(snd_f[label_combs[fig_id][0]][unit_name]['p_value'])
                    label1 += f" ({MRL:.2f}; {pval2text(pv)})"
                if label_combs[fig_id][1] in snd_f:
                    MRL = np.array(snd_f[label_combs[fig_id][1]][unit_name]['MRL_real'])
                    pv  = np.array(snd_f[label_combs[fig_id][1]][unit_name]['p_value'])
                    label2 += f" ({MRL:.2f}; {pval2text(pv)})"

        ax = fig.add_subplot(rows, cols, i+1)
        ax.hist(bins[:-1], bins=bins, weights=psth1, edgecolor='black', color=color_combs[fig_id][0], alpha=0.7, label=label1)
        ax.hist(bins[:-1], bins=bins, weights=psth2, edgecolor='black', color=color_combs[fig_id][1], alpha=0.7, label=label2)
        ax.axvline(0, color='black', ls='--')
        ax.axvspan(0, bgr_dur, alpha=0.3, color='gray')
        ax.axvspan(0 - hw, 0 - hw + bgr_dur, alpha=0.3, color='gray')
        ax.set_title(unit_name, fontsize=14)
        ax.legend(loc='lower right', prop={'size': 10})
        ax.set_xlim(-hw, hw)
        if i % 3 == 0:
            ax.set_ylabel("Firing Rate, Hz", fontsize=14)
            
    fig.tight_layout()
    fig.savefig(snakemake.output[fig_id])


# distractor figures
if (len(idxs_di1_ev) < 5 or len(idxs_di2_ev) < 5) and ((not (len(idxs_di1_ev) == 0 or len(idxs_di2_ev) == 0)) or 0 < len(idxs_di1_ev) < 5 and 0 < len(idxs_di2_ev) < 5):
    # not enough samples, write empty file
    f_name = os.path.join(os.path.dirname(snakemake.output[0]), 'psth_distractors.pdf')
    fig = plt.figure(figsize=(4, 4))
    fig.savefig(f_name)

else:
    label_combs = [
        ['TGT', 'Dis1'],
        ['TGT', 'Dis2'],
        ['Dis1', 'Dis2']
    ]

    color_combs = [
        ['tab:orange', 'navy'],
        ['tab:orange', 'green'],
        ['navy', 'green'],
    ]

    # saving in batches of 150 otherwise image is too large
    batch_size = 100
    batch_count = int(np.ceil(len(units_to_plot)/batch_size))

    for k in range(batch_count):
        units_selected = units_to_plot[k*batch_size:(k+1)*batch_size]

        rows = len(units_selected)
        cols = 3
        fig = plt.figure(figsize=(4*cols, 4*rows))

        for i, unit_name in enumerate(units_selected):
            for j, (idxs_ev_1, idxs_ev_2) in enumerate([(idxs_tgt_ev, idxs_di1_ev), (idxs_tgt_ev, idxs_di2_ev), (idxs_di1_ev, idxs_di2_ev)]):

                bins, psth1 = get_spike_counts(spike_times[unit_name], sound_events[idxs_ev_1][:, 0], hw=hw, bin_count=bc)
                bins, psth2 = get_spike_counts(spike_times[unit_name], sound_events[idxs_ev_2][:, 0], hw=hw, bin_count=bc)
                
                ax = fig.add_subplot(rows, cols, 3*i + j + 1)
                ax.hist(bins[:-1], bins=bins, weights=psth1, edgecolor='black', color=color_combs[j][0], alpha=0.7, label=label_combs[j][0])
                ax.hist(bins[:-1], bins=bins, weights=psth2, edgecolor='black', color=color_combs[j][1], alpha=0.7, label=label_combs[j][1])
                ax.axvline(0, color='black', ls='--')
                ax.axvspan(0, bgr_dur, alpha=0.3, color='gray')
                ax.axvspan(0 - hw, 0 - hw + bgr_dur, alpha=0.3, color='gray')
                ax.set_title(unit_name, fontsize=14)
                ax.legend(loc='lower right', prop={'size': 10})
                ax.set_xlim(-hw, hw)
                if i % 3 == 0:
                    ax.set_ylabel("Firing Rate, Hz", fontsize=14)
                
        f_name = os.path.join(os.path.dirname(snakemake.output[0]), 'psth_distractors_%s.pdf' % str(k+1))
        fig.tight_layout()
        fig.savefig(f_name)

if passive:
    #stim_combs = [(1, 2), (0, 1)]  # stimulus combinations to plot
    cols = 3
    rows = int(np.ceil(len(units_to_plot)/cols))

    hw = snakemake.config['psth']['micro']['latency']
    bc = snakemake.config['psth']['micro']['bin_count']

    speed_max = 0.04
    speed_ev = tl[sound_events[:, 2].astype(np.int32)][:, 3]
    idxs_sta_ev = np.where(speed_ev < speed_max)[0]
    idxs_run_ev = np.where(speed_ev > speed_max)[0]

    fig_titles = [
        "frequencies",
        "durations"
    ]

    stim_comb_idxs = [
        [idxs_ev for key, idxs_ev in idxs_ev_passive.items() if 'F' in key],  # frequencies plot
        [idxs_ev for key, idxs_ev in idxs_ev_passive.items() if 'D' in key]  # durations plot
        # ignore stationary / run for passive
        # [np.intersect1d(idxs_bgr_ev, idxs_sta_ev), np.intersect1d(idxs_tgt_ev, idxs_sta_ev)],  # BGR stationary / TGT stationary
        # [np.intersect1d(idxs_bgr_ev, idxs_sta_ev), np.intersect1d(idxs_bgr_ev, idxs_run_ev)],  # BGR stationary / BGR run
        # [np.intersect1d(idxs_sil_ev, idxs_sta_ev), np.intersect1d(idxs_sil_ev, idxs_run_ev)],  # SIL stationary / SIL run
    ]

    label_combs = [
        [f"{sound["freq"]}Hz" for key, sound in cfg['sound']['sounds'].items() if 'F' in key],
        [f"{float(sound["duration"])*1000}ms" for key, sound in cfg['sound']['sounds'].items() if 'D' in key]
        # ignore stationary / run for passive
        # ['bgr_sta', 'tgt_sta'],
        # ['bgr_sta', 'bgr_run'],
        # ['sil_sta', 'sil_run'],
    ]

    # colors for the different frequencies should be shades of red
    # colors for the different durations should be shades of blue
    # color_combs = [
    #     # colors for the different frequencies should be different shades of red, not the same color
    #     [plt.cm.Reds(i / max(1, (len([key for key in cfg['sound']['sounds'] if 'F' in key]) - 1))) 
    #         for i, key in enumerate([key for key in cfg['sound']['sounds'] if 'F' in key])],

    #     [plt.cm.Blues(i / max(1, (len([key for key in cfg['sound']['sounds'] if 'D' in key]) - 1))) 
    #         for i, key in enumerate([key for key in cfg['sound']['sounds'] if 'D' in key])]
    #     # ignore stationary / run for passive
    #     # ['tab:blue', 'tab:orange'],
    #     # ['navy', 'tab:blue'],
    #     # ['grey', 'tab:red'],
    # ]
    # different colors
    color_combs = [
        ['#d0d1e6','#a6bddb','#74a9cf','#3690c0','#0570b0','#034e7b'],
        ['#addd8e','#78c679','#41ab5d','#238443','#005a32'],
        # [f'C{i}' for i, key in enumerate([key for key in cfg['sound']['sounds'] if 'F' in key])],
        # [f'C{i}' for i, key in enumerate([key for key in cfg['sound']['sounds'] if 'D' in key])]
        # ignore stationary / run for passive
        # ['tab:blue', 'tab:orange'],
        # ['navy', 'tab:blue'],
        # ['grey', 'tab:red'],
    ]

    # TGT, BGR, SIL bar plot figures
    for fig_id, stim_comb in enumerate(stim_comb_idxs):

        # figure / file for each stimulus combination
        fig = plt.figure(figsize=(4*cols, 4*rows))

        for i, unit_name in enumerate(units_to_plot):
            ax = fig.add_subplot(rows, cols, i+1) # one subplot for each unit
            for ev_id, idxs_ev in enumerate(stim_comb):
                bins, psth = get_spike_counts(spike_times[unit_name], sound_events[idxs_ev][:, 0], hw=hw, bin_count=bc)

                # if sound phase locking exists - plot with the label
                label = label_combs[fig_id][ev_id]
                # TODO: ignore this for the moment
                # if os.path.exists(sound_phase_lock_file):
                #     with h5py.File(sound_phase_lock_file, 'r') as snd_f:
                #         if label_combs[fig_id][0] in snd_f:
                #             MRL = np.array(snd_f[label_combs[fig_id][0]][unit_name]['MRL_real'])
                #             pv  = np.array(snd_f[label_combs[fig_id][0]][unit_name]['p_value'])
                #             label1 += f" ({MRL:.2f}; {pval2text(pv)})"
                #         if label_combs[fig_id][1] in snd_f:
                #             MRL = np.array(snd_f[label_combs[fig_id][1]][unit_name]['MRL_real'])
                #             pv  = np.array(snd_f[label_combs[fig_id][1]][unit_name]['p_value'])
                #             label2 += f" ({MRL:.2f}; {pval2text(pv)})"

                print(ev_id, label, color_combs[fig_id][ev_id])
                ax.hist(bins[:-1], bins=bins, weights=psth, color=color_combs[fig_id][ev_id], histtype='step', label=label, linewidth=3)


            ax.axvline(0, color='black', ls='--')
            ax.axvspan(0, bgr_dur, alpha=0.3, color='gray')
            ax.axvspan(0 - hw, 0 - hw + bgr_dur, alpha=0.3, color='gray')
            ax.set_title(f"{unit_name} ({depth[unit_name]}um)", fontsize=14)
            ax.legend(loc='lower right', prop={'size': 10})
            # ax.set_xlim(-hw, hw)
            ax.set_xlim(0, hw)
            if i % 3 == 0:
                ax.set_ylabel("Firing Rate, Hz", fontsize=14)
                
        fig.tight_layout()
        # TODO: this is done in an ugly way
        f_name = os.path.join(os.path.dirname(snakemake.output[0]), f'{fig_titles[fig_id]}.pdf')
        fig.savefig(f_name)
