import os, sys
import h5py
import json
import numpy as np

# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.decoders import run_decoder_with_shuffles

aep_type    = snakemake.config['state_decoder']['aep_type']
kernel_size = snakemake.config['state_decoder']['kernel_size']

meta_file = snakemake.input[0]
segm_file = snakemake.input[1]
nmap_file = snakemake.input[2]
aepm_file = snakemake.input[3]

# ------ reading LFP / population inputs ------

with h5py.File(meta_file, 'r') as f:
    tl = np.array(f['processed']['timeline'])
    sound_events = np.array(f['processed']['sound_events'])
    tgt_mx = np.array(f['processed']['target_matrix'])

with h5py.File(nmap_file, 'r') as f:
    pop_ev = np.array(f['response_manifold'])[:, 0]
    pop_su = np.array(f['response_manifold'])[:, 1]

with h5py.File(aepm_file, 'r') as f:
    lfp_ev = np.array(f[aep_type][kernel_size])[:, 0]
    lfp_su = np.array(f[aep_type][kernel_size])[:, 1]

# reading state indices
state_ids    = ['idxs_tgt_sta_succ', 'idxs_bgr_sta', 'idxs_bgr_run', 'idxs_sil_sta', 'idxs_sil_run']
state_labels = ['Target', 'BGR/sta', 'BGR/run', 'No stim./sta', 'No stim./run']
state_colors = ['tab:orange', 'tab:blue', 'red', 'grey', 'black']
with h5py.File(segm_file, 'r') as f:
    state_idxs = {}
    for idxs_name in state_ids:
        state_idxs[idxs_name] = np.array(f[idxs_name]).astype(np.int32)

# ------ building state labels / feature vectors ------

# 0 - Target, 1 - Background STA, 2 - Background RUN, 3 - No stimulus STA, 4 - No stimulus RUN

labels = -1 * np.ones(len(sound_events))
labels[state_idxs['idxs_tgt_sta_succ']] = 0  # Target
labels[state_idxs['idxs_bgr_sta']] = 1        # Background STA
labels[state_idxs['idxs_bgr_run']] = 2        # Background RUN
labels[state_idxs['idxs_sil_sta']] = 3        # No stimulus - STA
labels[state_idxs['idxs_sil_run']] = 4        # No stimulus RUN

# filter unknown / boundary states
known = labels >= 0
lfp_ev = lfp_ev[known]
lfp_su = lfp_su[known]
pop_ev = pop_ev[known]
pop_su = pop_su[known]

features = {
    'X_lfp_ev':       lfp_ev[:, None],
    'X_lfp_su':       lfp_su[:, None],
    'X_pop_ev':       pop_ev[:, None],
    'X_pop_su':       pop_su[:, None],
    'X_lfp_comb':     np.column_stack([lfp_ev, lfp_su]),
    'X_pop_comb':     np.column_stack([pop_ev, pop_su]),
    'X_lfp_pop_comb': np.column_stack([lfp_ev, lfp_su, pop_ev, pop_su])
}

# ------ do decoding with shuffles -------

results = {}
for name, X in features.items():
    (
        acc_real,
        cm_real,
        y_real,
        y_pred_real,
        acc_shuf,
        cm_shuf_mean,
        cm_shuf_std,
    ) = run_decoder_with_shuffles(X, labels[known])

    with h5py.File(snakemake.output[0], 'a') as f:
        if name in f:
            del f[name]

        grp = f.create_group(name)
        grp.create_dataset('acc_real', data=acc_real)
        grp.create_dataset('cm_real', data=cm_real)
        grp.create_dataset('y_real', data=y_real)
        grp.create_dataset('y_pred_real', data=y_pred_real)
        grp.create_dataset('acc_shuf', data=acc_shuf)
        grp.create_dataset('cm_shuf_mean', data=cm_shuf_mean)
        grp.create_dataset('cm_shuf_std', data=cm_shuf_std)
