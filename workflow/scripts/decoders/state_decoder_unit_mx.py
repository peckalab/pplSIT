import os, sys
import h5py
import json
import numpy as np

# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.decoders import run_decoder_with_shuffles

electrodes = snakemake.config['state_decoder']['unit_mx_electrodes']  # filter electrodes if needed

meta_file = snakemake.input[0]
segm_file = snakemake.input[1]
actm_file = snakemake.input[2]

# ------ reading LFP / population inputs ------

with h5py.File(meta_file, 'r') as f:
    tl = np.array(f['processed']['timeline'])
    sound_events = np.array(f['processed']['sound_events'])
    tgt_mx = np.array(f['processed']['target_matrix'])

with h5py.File(actm_file, 'r') as f:
    unit_mx = np.array(f['mx_250ms']['mx'])  # units x time bins (1 value less than sound events)
unit_mx = np.hstack([unit_mx, unit_mx[:, -1][:, None]]).T  # duplicate last values and transpose (time bins x units)

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
X = unit_mx[known]
y = labels[known]

# ------ do decoding with shuffles -------

(
    acc_real,
    cm_real,
    y_real,
    y_pred_real,
    acc_shuf,
    cm_shuf_mean,
    cm_shuf_std,
) = run_decoder_with_shuffles(X, y)

name = 'X_unit_mx'
with h5py.File(snakemake.output[0], 'w') as f:
    grp = f.create_group(name)
    grp.create_dataset('acc_real', data=acc_real)
    grp.create_dataset('cm_real', data=cm_real)
    grp.create_dataset('y_real', data=y_real)
    grp.create_dataset('y_pred_real', data=y_pred_real)
    grp.create_dataset('acc_shuf', data=acc_shuf)
    grp.create_dataset('cm_shuf_mean', data=cm_shuf_mean)
    grp.create_dataset('cm_shuf_std', data=cm_shuf_std)
