import h5py, os, sys, json
import numpy as np
import itertools

# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.states import get_state_as_periods


# read datasets
meta_file = snakemake.input[0]
unit_file = snakemake.input[1]
actm_file = snakemake.input[2]
segm_file = snakemake.input[3]
smk_cfg   = snakemake.config['vectors']

s_path  = os.path.dirname(meta_file)
session = os.path.basename(s_path)
animal  = session.split('_')[0]

with h5py.File(meta_file, 'r') as f:
    tl = np.array(f['processed']['timeline'])
    tgt_mx = np.array(f['processed']['target_matrix'])
    sound_events = np.array(f['processed']['sound_events'])
    cfg = json.loads(f['processed'].attrs['parameters'])

unit_mfrs = []
with h5py.File(unit_file, 'r') as f:
    unit_ids = ([name for name in f])
    for unit_name in unit_ids:
        unit_mfrs.append(float(np.array(f[unit_name]['mean_firing_rate'])))
unit_mfrs = np.array(unit_mfrs)

# getting experimental states
with h5py.File(segm_file, 'r') as f:
    tgt_sta_succ_mx = np.array(f['tgt_sta_succ_mx'])
    tgt_sta_mx = np.array(f['tgt_sta_mx'])
    bgr_sta_mx = np.array(f['bgr_sta_mx'])
    sil_sta_mx = np.array(f['sil_sta_mx'])
tgt_sta_mx_t3 = tgt_sta_mx.copy()
tgt_sta_mx_t3[:, 0] = tgt_sta_mx_t3[:, 1] - 12  # last 3 seconds in target

# filtering units
if session in smk_cfg['units_to_exclude']:
    idxs_special = [unit_ids.index(x) for x in smk_cfg['units_to_exclude'][session]]
else:
    idxs_special = []
idxs_mfr_min = np.where(unit_mfrs > smk_cfg['min_mfr'])[0]

idxs_units_filt = np.setdiff1d(idxs_mfr_min, idxs_special)

# unit activity matrix
with h5py.File(actm_file, 'r') as f:
    u_group = smk_cfg['unit_mx']['group']
    u_type  = smk_cfg['unit_mx']['type']
    u_step  = smk_cfg['unit_mx']['step']
    unit_mx_ev = np.array(f[u_group][u_type])[:, ::u_step].T
unit_mx = unit_mx_ev.T[idxs_units_filt].T  # filter required units



results = {}

# ----- correlations WITHIN the same state -----
state_names = ['tgt_sta', 'bgr_sta', 'sil_sta']

# all matrices are periods of: first pulse, last pulse, pos X, pos Y
state_mxs  = [tgt_sta_mx, bgr_sta_mx, sil_sta_mx]
pulses_rel = [
    [x for x in range(-12, 40)],
    [x for x in range(-12, 12)],
    [x for x in range(-12, 12)]
]

for i, state_mx in enumerate(state_mxs):
    corr_medians = []
    corr_tgt_mx = []
    for rel_pulse_idx in pulses_rel[i]:
        combs = list(itertools.combinations(range(len(state_mx)), 2))  # combinations of vector IDs

        for idx_1, idx_2 in combs:
            idx_p1 = int(state_mx[idx_1][0]) + rel_pulse_idx
            idx_p2 = int(state_mx[idx_2][0]) + rel_pulse_idx
            if idx_p1 >= len(unit_mx) or idx_p2 >= len(unit_mx):
                continue
            
            v1 = unit_mx[int(state_mx[idx_1][0]) + rel_pulse_idx]
            v2 = unit_mx[int(state_mx[idx_2][0]) + rel_pulse_idx]

            corr = np.corrcoef(v1, v2)[0][1]
            dist = np.sqrt( (state_mx[idx_1][2] - state_mx[idx_2][2])**2 + (state_mx[idx_1][3] - state_mx[idx_2][3])**2 )
            time = np.abs(state_mx[idx_2][0] - state_mx[idx_1][0]) / 4  # in seconds
                
            # pulse id rel to target onset, corr coeff, dist b/w targets, time b/w targets
            corr_tgt_mx.append(np.array([rel_pulse_idx, corr, dist, time]))
    results[state_names[i]] = np.array(corr_tgt_mx)


# ----- correlations BETWEEN states -----
state_names = ['bgr_sta-tgt_sta', 'bgr_sta-tgt_sta_t3']

state_combs = [
    [bgr_sta_mx, tgt_sta_mx],
    [bgr_sta_mx, tgt_sta_mx_t3],
]
pulses_rel = [
    [x for x in range(-12, 12)],
    [x for x in range(-12, 12)]
]

for i, (sel_mx1, sel_mx2) in enumerate(state_combs):
    corr_medians = []
    corr_tgt_mx = []
    for rel_pulse_idx in pulses_rel[i]:
        combs = list(itertools.product(np.arange(len(sel_mx1)), np.arange(len(sel_mx2))))  # combinations of vector IDs

        for idx_1, idx_2 in combs:
            v1 = unit_mx[int(sel_mx1[idx_1][0]) + rel_pulse_idx]
            v2 = unit_mx[int(sel_mx2[idx_2][0]) + rel_pulse_idx]

            corr = np.corrcoef(v1, v2)[0][1]
            dist = np.sqrt( (sel_mx1[idx_1][2] - sel_mx2[idx_2][2])**2 + (sel_mx1[idx_1][3] - sel_mx2[idx_2][3])**2 )
            time = np.abs(sel_mx2[idx_2][0] - sel_mx1[idx_1][0]) / 4  # in seconds
                
            # pulse id rel to target onset, corr coeff, dist b/w targets, time b/w targets
            corr_tgt_mx.append(np.array([rel_pulse_idx, corr, dist, time]))
    results[state_names[i]] = np.array(corr_tgt_mx)

# dump to H5
with h5py.File(snakemake.output[0], 'w') as out_file:
    for name, corr_mx in results.items():
        out_file.create_dataset(name, data=corr_mx)