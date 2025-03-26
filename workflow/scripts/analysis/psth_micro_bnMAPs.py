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
# filter units for publications, if needed
special_units = snakemake.config['PSTH_micro_bnMAPs']['units']
if len(special_units) > 0:
    units_to_plot = [unit for unit in units_to_plot if unit in special_units]


with h5py.File(snakemake.input[2], 'r') as f:
    speed = np.array(f['speed'])
    hd    = np.array(f['hd'])

# with h5py.File(snakemake.input[3], 'r') as f:
#     nfit = np.array(f['response_manifold'])
#     unit_mx_ev = np.array(f['unit_mx_proc_ev'])
#     unit_mx_su = np.array(f['unit_mx_proc_su'])

with h5py.File(snakemake.input[3], 'r') as f:
    idxs_tgt_succ_state_ev = np.array(f['idxs_tgt_succ_state_ev'])


# configuration
event_types = [0, 1, 2, -1]  # SIL, BGR, TGT, NOI - order matters
colors = {0: 'gray', 1: 'tab:blue', 2: 'tab:orange', -1: 'red'}
ev_names = {0: 'SIL', 1: 'BGR', 2: 'TGT', -1: 'NOI'}
bgr_dur    = cfg['sound']['sounds']['background']['duration']  # in seconds
#stim_combs = [(1, 2), (0, 1)]  # stimulus combinations to plot
cols = 3
rows = int(np.ceil(len(units_to_plot)/cols))
colors = {0: 'gray', 1: 'tab:blue', 2: 'tab:orange', -1: 'red'}
ev_names = {0: 'SIL', 1: 'BGR', 2: 'TGT', -1: 'NOI'}
#latency = cfg['sound']['latency']

hw = snakemake.config['psth']['micro']['latency']
bc = snakemake.config['psth']['micro']['bin_count']
figsize = snakemake.config['psth']['micro']['figsize']

speed_max = 0.04
#speed_ev = tl[sound_events[:, 2].astype(np.int32)][:, 3]  # speed taken from timeline
speed_ev = speed[sound_events[:, 2].astype(np.int32)]  # better - speed taken from DLC
idxs_sta_ev = np.where(speed_ev < speed_max)[0]
idxs_run_ev = np.where(speed_ev > speed_max)[0]
idxs_bgr_ev = np.where(sound_events[:, 1] == 1)[0]
idxs_sil_ev = np.where(sound_events[:, 1] == 0)[0]
idxs_tgt_ev = np.where(sound_events[:, 1] == 2)[0]
idxs_di1_ev = np.where(sound_events[:, 1] == 3)[0]
idxs_di2_ev = np.where(sound_events[:, 1] == 4)[0]  # for the moment limit to 2 distractors

idxs_bU_ev = np.array([x for x in range(len(sound_events)) if not x in idxs_tgt_succ_state_ev])

# multiplexing events: sound, locomotion, engagement

# ------ NO STIMULUS ------

# 000: in silence, stationary non-engaged condition
idxs_000 = np.intersect1d(idxs_sta_ev, idxs_sil_ev)
idxs_000 = np.intersect1d(idxs_000, idxs_bU_ev)

# 001: in silence, stationary engaged condition
idxs_001 = np.intersect1d(idxs_sta_ev, idxs_sil_ev)
idxs_001 = np.intersect1d(idxs_001, idxs_tgt_succ_state_ev)

# 010: in silence, running non-engaged condition
idxs_010 = np.intersect1d(idxs_run_ev, idxs_sil_ev)
idxs_010 = np.intersect1d(idxs_010, idxs_bU_ev)

# 011: in silence, running engaged - should not happen!
idxs_011 = np.intersect1d(idxs_run_ev, idxs_sil_ev)
idxs_011 = np.intersect1d(idxs_011, idxs_tgt_succ_state_ev)

# ------ BGR STIMULUS ------

# 100: in trial, stationary non-engaged condition
idxs_100 = np.intersect1d(idxs_sta_ev, idxs_bgr_ev)
idxs_100 = np.intersect1d(idxs_100, idxs_bU_ev)

# 101: in trial, stationary engaged condition
idxs_101 = np.intersect1d(idxs_sta_ev, idxs_bgr_ev)
idxs_101 = np.intersect1d(idxs_101, idxs_tgt_succ_state_ev)

# 110: in trial, running non-engaged condition
idxs_110 = np.intersect1d(idxs_run_ev, idxs_bgr_ev)
idxs_110 = np.intersect1d(idxs_110, idxs_bU_ev)

# 111: in trial, running engaged - should not happen!
idxs_111 = np.intersect1d(idxs_run_ev, idxs_bgr_ev)
idxs_111 = np.intersect1d(idxs_111, idxs_tgt_succ_state_ev)

# ------ TGT STIMULUS ------

# 200: in TGT, stationary non-engaged condition
idxs_200 = np.intersect1d(idxs_sta_ev, idxs_tgt_ev)
idxs_200 = np.intersect1d(idxs_200, idxs_bU_ev)

# 201: in TGT, stationary engaged condition
idxs_201 = np.intersect1d(idxs_sta_ev, idxs_tgt_ev)
idxs_201 = np.intersect1d(idxs_201, idxs_tgt_succ_state_ev)

# 210: in TGT, running non-engaged condition
idxs_210 = np.intersect1d(idxs_run_ev, idxs_tgt_ev)
idxs_210 = np.intersect1d(idxs_210, idxs_bU_ev)

# 211: in TGT, running engaged - should not happen!
idxs_211 = np.intersect1d(idxs_run_ev, idxs_tgt_ev)
idxs_211 = np.intersect1d(idxs_211, idxs_tgt_succ_state_ev)


stim_comb_idxs = [
    [idxs_000, idxs_001],  # SIL sta bU / SIL sta bE
    [idxs_010, idxs_000],  # SIL run bU / SIL sta bU
    [idxs_100, idxs_101],  # BGR sta bU / BGR sta bE
    [idxs_110, idxs_101],  # BGR run bU / BGR sta bU
]

label_combs = [
    ['Stationary', 'Stationary, Engaged'],
    ['Running', 'Stationary'],
    ['Stationary', 'Stationary, Engaged'],
    ['Running', 'Stationary']
]

color_combs = [
    ['grey', 'tab:green'],
    ['tab:red', 'grey'],
    ['tab:blue', 'tab:green'],
    ['tab:red', 'tab:blue']
]

# TGT, BGR, SIL bar plot figures
for fig_id, stim_comb in enumerate(stim_comb_idxs):
    idxs_ev_1 = stim_comb[0]
    idxs_ev_2 = stim_comb[1]

    # figure / file for each stimulus combination
    fig = plt.figure(figsize=(figsize*cols, figsize*rows))

    for i, unit_name in enumerate(units_to_plot):
        bins, psth1 = get_spike_counts(spike_times[unit_name], sound_events[idxs_ev_1][:, 0], hw=hw, bin_count=bc)
        bins, psth2 = get_spike_counts(spike_times[unit_name], sound_events[idxs_ev_2][:, 0], hw=hw, bin_count=bc)
        
        ax = fig.add_subplot(rows, cols, i+1)
        ax.hist(bins[:-1], bins=bins, weights=psth1, edgecolor='black', color=color_combs[fig_id][0], alpha=0.9, label=label_combs[fig_id][0])
        ax.hist(bins[:-1], bins=bins, weights=psth2, edgecolor='black', color=color_combs[fig_id][1], alpha=0.9, label=label_combs[fig_id][1])
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