import h5py, os, sys, json
import numpy as np

# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.CCR import xcorr_session_with_shuffle

#state_labels = ['Target', 'Background (sta)', 'Background (run)', 'No stimulus (sta)', 'No stimulus (run)']
#state_colors = ['tab:orange', 'mediumblue', 'deepskyblue', 'black', 'grey']
state_ids    = ['idxs_tgt_sta_succ', 'idxs_bgr_sta', 'idxs_bgr_run', 'idxs_sil_sta', 'idxs_sil_run']
aep_type = 'avg_across_channels'
kernel_size = '4'

s_path  = os.path.dirname(snakemake.input[0])
session = os.path.basename(s_path)
animal  = session.split('_')[0]

idxs_states = {}
with h5py.File(snakemake.input[1], 'r') as f:
    for st_name in state_ids:
        idxs_states[st_name] = np.array(f[st_name])
with h5py.File(snakemake.input[2], 'r') as f:
    EV_aep = np.array(f[aep_type][kernel_size][:, 0])
    SU_aep = np.array(f[aep_type][kernel_size][:, 1])
with h5py.File(snakemake.input[3], 'r') as f:
    response_manifold = np.array(f['response_manifold'])
    EV_pop = response_manifold[:, 0]
    SU_pop = response_manifold[:, 1]

#res = xcorr_session(m1, m2, cond_indices=idxs_states, fs_hz=4, max_lag_s=20)
res_aep = xcorr_session_with_shuffle(EV_aep, SU_aep, idxs_states, fs_hz=4, max_lag_s=20, n_shuffle=100, shuffle_mode='cshift', seed=42)
res_pop = xcorr_session_with_shuffle(EV_pop, SU_pop, idxs_states, fs_hz=4, max_lag_s=20, n_shuffle=100, shuffle_mode='cshift', seed=42)

# dump results
with h5py.File(snakemake.output[0], 'w') as f:
    LFP_grp = f.create_group('LFP')
    pop_grp = f.create_group('population')

    for cond_id, cond_res in res_aep.items():
        grp = LFP_grp.create_group(cond_id)
        for res_key, data in cond_res.items():
            grp.create_dataset(res_key, data=data)

    for cond_id, cond_res in res_pop.items():
        grp = pop_grp.create_group(cond_id)
        for res_key, data in cond_res.items():
            grp.create_dataset(res_key, data=data)

