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


# reading some configs
with h5py.File(snakemake.input[0], 'r') as f:
    cfg = json.loads(f['processed'].attrs['parameters'])
    sound_events = np.array(f['processed']['sound_events'])
    tl = np.array(f['processed']['timeline'])

spike_times = {}
with h5py.File(snakemake.input[1], 'r') as f:
    units_to_plot = get_unit_names_sorted([name for name in f])
    for unit_name in units_to_plot:
        spike_times[unit_name] = np.array(f[unit_name]['spike_times'])

with h5py.File(snakemake.input[2], 'r') as f:
    ensemble_tl = np.array(f['AL'])
ensemble_ev = ensemble_tl[sound_events[:, 2].astype(np.int32)]

# configuration
event_types = [0, 1, 2, -1]  # SIL, BGR, TGT, NOI - order matters
colors = {0: 'gray', 1: 'tab:blue', 2: 'tab:orange', -1: 'red'}
ev_names = {0: 'SIL', 1: 'BGR', 2: 'TGT', -1: 'NOI'}
bgr_dur    = cfg['sound']['sounds']['background']['duration']  # in seconds
cols = 3
rows = int(np.ceil(len(units_to_plot)/cols))
colors = {0: 'gray', 1: 'tab:blue', 2: 'tab:orange', -1: 'red'}
ev_names = {0: 'SIL', 1: 'BGR', 2: 'TGT', -1: 'NOI'}

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

idxs_AL_ev  = np.where(ensemble_ev > 0)[0]
idxs_PH_ev  = np.where(ensemble_ev < 0)[0]

# states from neuronal ensembles
# - BGR sta AL / BGR sta PH
# - BGR sta AL / TGT
idxs_bgr_sta_ev = np.intersect1d(idxs_bgr_ev, idxs_sta_ev)
idxs_bgr_run_ev = np.intersect1d(idxs_bgr_ev, idxs_sta_ev)
idxs_bgr_sta_AL_ev = np.intersect1d(idxs_bgr_sta_ev, idxs_AL_ev)
idxs_bgr_sta_PH_ev = np.intersect1d(idxs_bgr_sta_ev, idxs_PH_ev)
idxs_tgt_AL_ev = np.intersect1d(idxs_tgt_ev, idxs_AL_ev)
idxs_tgt_PH_ev = np.intersect1d(idxs_tgt_ev, idxs_PH_ev)

stim_comb_idxs = [
    [idxs_bgr_sta_PH_ev, idxs_bgr_sta_AL_ev],  # BGR stationary AL / PH
    [idxs_tgt_ev, idxs_bgr_sta_AL_ev],  # BGR stationary AL / TGT
    [idxs_tgt_PH_ev, idxs_tgt_AL_ev],  # TGT AL / PH
]

label_combs = [
    ['BGR sta PH', 'BGR sta AL'],
    ['TGT', 'BGR sta AL'],
    ['TGT PH', 'TGT AL'],
]

color_combs = [
    ['grey', 'navy'],
    ['tab:orange', 'navy'],
    ['tab:orange', 'red'],
]

# TGT, BGR, SIL bar plot figures
for fig_id, stim_comb in enumerate(stim_comb_idxs):
    idxs_ev_1 = stim_comb[0]
    idxs_ev_2 = stim_comb[1]
    label1 = "%s (%d)" % (label_combs[fig_id][0], len(idxs_ev_1))
    label2 = "%s (%d)" % (label_combs[fig_id][1], len(idxs_ev_2))

    # figure / file for each stimulus combination
    fig = plt.figure(figsize=(4*cols, 4*rows))

    for i, unit_name in enumerate(units_to_plot):
        bins, psth1 = get_spike_counts(spike_times[unit_name], sound_events[idxs_ev_1][:, 0], hw=hw, bin_count=bc)
        bins, psth2 = get_spike_counts(spike_times[unit_name], sound_events[idxs_ev_2][:, 0], hw=hw, bin_count=bc)
        label1 = "%s (%d)" % (label_combs[fig_id][0], len(idxs_ev_1))
        label2 = "%s (%d)" % (label_combs[fig_id][1], len(idxs_ev_2))

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
