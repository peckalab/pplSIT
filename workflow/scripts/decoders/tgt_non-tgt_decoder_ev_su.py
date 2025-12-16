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
nmap_file = snakemake.input[1]
aepm_file = snakemake.input[2]

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

# ------ building state labels / feature vectors ------

idxs_bgr = np.where(sound_events[:, 1] == 1)[0]
idxs_tgt = np.where(sound_events[:, 1] == 2)[0]
idxs_sound_on = np.union1d(idxs_bgr, idxs_tgt)
idxs_all = np.arange(len(sound_events))

labels = np.zeros(len(sound_events))
labels[idxs_tgt] = 1  # Target vs. non-Target

features = {
    'X_lfp_ev':       lfp_ev[:, None],
    'X_lfp_su':       lfp_su[:, None],
    'X_pop_ev':       pop_ev[:, None],
    'X_pop_su':       pop_su[:, None],
    'X_ev_comb':      np.column_stack([lfp_ev, pop_ev]),
    'X_su_comb':      np.column_stack([lfp_su, pop_su]),
    'X_lfp_comb':     np.column_stack([lfp_ev, lfp_su]),
    'X_pop_comb':     np.column_stack([pop_ev, pop_su]),
}

groups = {
    'all': idxs_all,
    'sound_on': idxs_sound_on,
}

# ------ do decoding with shuffles -------

results = {}
for grp_name, idxs in groups.items():
    with h5py.File(snakemake.output[0], 'a') as f:
        if grp_name in f:
            del f[name]
        grp_f = f.create_group(grp_name)
            
    for name, X in features.items():
        # decoder with shuffles
        (
            acc_real,
            cm_real,
            y_real,
            y_pred_real,
            acc_shuf,
            cm_shuf_mean,
            cm_shuf_std,
        ) = run_decoder_with_shuffles(X[idxs], labels[idxs])

        with h5py.File(snakemake.output[0], 'a') as f:
            grp_f = f[grp_name]
            if name in grp_f:
                del grp_f[name]

            grp = grp_f.create_group(name)
            grp.create_dataset('acc_real', data=acc_real)
            grp.create_dataset('cm_real', data=cm_real)
            grp.create_dataset('y_real', data=y_real)
            grp.create_dataset('y_pred_real', data=y_pred_real)
            grp.create_dataset('acc_shuf', data=acc_shuf)
            grp.create_dataset('cm_shuf_mean', data=cm_shuf_mean)
            grp.create_dataset('cm_shuf_std', data=cm_shuf_std)
