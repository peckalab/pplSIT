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

# getting standard experimental states
mx_types = ['tgt_sta_mx', 'tgt_sta_succ_mx', 'bgr_sta_mx', 'sil_sta_mx', 'bgr_run_mx', 'sil_run_mx', 'bgr_vis_mx', 'sil_vis_mx']

with h5py.File(segm_file, 'r') as f:
    event_matrices = {}
    for mx_type in mx_types:
        event_matrices[mx_type] = np.array(f[mx_type])

# adding last 3 seconds in scores
tgt_sta_mx_t3 = event_matrices['tgt_sta_mx'].copy()
tgt_sta_mx_t3[:, 0] = tgt_sta_mx_t3[:, 1] - 12  # last 3 seconds in target
mx_types.append('tgt_sta_mx_t3')
event_matrices['tgt_sta_mx_t3'] = tgt_sta_mx_t3

# adding AL / PH states independent on sound
with h5py.File(segm_file, 'r') as f:
    if 'bgr_sta_AL_mx' in f:
        bgr_sta_AL_mx = np.array(f['bgr_sta_AL_mx'])
        sil_sta_AL_mx = np.array(f['sil_sta_AL_mx'])
        bgr_sta_PH_mx = np.array(f['bgr_sta_PH_mx'])
        sil_sta_PH_mx = np.array(f['sil_sta_PH_mx'])
        mx_types.append('AL')
        mx_types.append('PH')
        event_matrices['AL'] = np.vstack([bgr_sta_AL_mx, sil_sta_AL_mx])
        event_matrices['PH'] = np.vstack([bgr_sta_PH_mx, sil_sta_PH_mx])

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

# taking mean pulse vector for visits, stays and last 3 seconds for target
event_vectors = {}
for mx_type in mx_types:
    vectors = []
    for event in event_matrices[mx_type]:
        act_mx_slice = unit_mx[int(event[0]):int(event[1])]
        vectors.append(act_mx_slice.mean(axis=0))
    event_vectors[mx_type] = np.array(vectors)


# ----- correlations WITHIN the same state -----

corr_mxs = {}
for ev_type in mx_types:
    vectors = event_vectors[ev_type]
    visits  = event_matrices[ev_type]
    combs = list(itertools.combinations(range(len(vectors)), 2))  # combinations of vector IDs

    corr_mx = []
    for idx_1, idx_2 in combs:
        v1 = vectors[idx_1]
        v2 = vectors[idx_2]

        corr = np.corrcoef(v1, v2)[0][1]
        dist = np.sqrt( (visits[idx_1][2] - visits[idx_2][2])**2 + (visits[idx_1][3] - visits[idx_2][3])**2 )
        time = np.abs(visits[idx_2][0] - visits[idx_1][0]) / 4  # in seconds

        corr_mx.append(np.array([corr, dist, time]))
    corr_mxs[ev_type] = np.array(corr_mx)

# ----- correlations BETWEEN states -----

event_combs = [
    ['bgr_vis_mx', 'tgt_sta_mx'],
    ['sil_vis_mx', 'tgt_sta_mx'],
    ['bgr_sta_mx', 'tgt_sta_mx'],
    ['sil_sta_mx', 'tgt_sta_mx'],
]

for event_comb in event_combs:
    vectors1 = event_vectors[event_comb[0]]
    vectors2 = event_vectors[event_comb[1]]
    visits1  = event_matrices[event_comb[0]]
    visits2  = event_matrices[event_comb[1]]
    combs = list(itertools.product(np.arange(len(vectors1)), np.arange(len(vectors2))))  # combinations of vector IDs

    corr_mx = []
    for idx_1, idx_2 in combs:
        v1 = vectors1[idx_1]
        v2 = vectors2[idx_2]

        corr = np.corrcoef(v1, v2)[0][1]
        dist = np.sqrt( (visits1[idx_1][2] - visits2[idx_2][2])**2 + (visits1[idx_1][3] - visits2[idx_2][3])**2 )
        time = np.abs(visits2[idx_2][0] - visits1[idx_1][0]) / 4  # in seconds

        corr_mx.append(np.array([corr, dist, time]))
    corr_mxs[f'{event_comb[0]}-{event_comb[1]}'] = np.array(corr_mx)

# dump to H5
with h5py.File(snakemake.output[0], 'w') as out_file:
    for name, corr_mx in corr_mxs.items():
        out_file.create_dataset(name, data=corr_mx)