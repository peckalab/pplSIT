import h5py, os, sys, json
import numpy as np
from scipy import stats, signal
from sklearn import decomposition
from sklearn.preprocessing import minmax_scale


# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.population import unit_activity_matrix
from utils.psth import get_psth_matrix, staple_spike_times
from utils.spiketrain import smooth_gaussian, smooth_rectangular


s_path  = os.path.dirname(snakemake.input[0])
session = os.path.basename(s_path)
animal  = session.split('_')[0]
electrodes = snakemake.config['nMAP_electrodes'][animal]  # electrodes in A1

smooth_su_size = snakemake.config['nMAP_EV_SU']['smooth_su_size']
smooth_ev_size = snakemake.config['nMAP_EV_SU']['smooth_ev_size']

# reading data
with h5py.File(snakemake.input[0], 'r') as f:
    sound_events = np.array(f['processed']['sound_events'])
    tl = np.array(f['processed']['timeline'])
    tgt_mx = np.array(f['processed']['target_matrix'])
    cfg = json.loads(f['processed'].attrs['parameters'])

with h5py.File(snakemake.input[1], 'r') as f:
    ev_bin_count = int(np.array(f['ev_bin_count']))
    ev_periods = np.array(f['ev_periods'])
    ev_bins = np.array(f['ev_bins'])
    ev_unit_mx = np.array(f['ev_unit_mx'])  # original activity matrix, 10 ms bins

    su_bin_count = int(np.array(f['su_bin_count']))
    su_periods = np.array(f['su_periods'])
    su_bins = np.array(f['su_bins'])
    su_unit_mx = np.array(f['su_unit_mx'])  # original activity matrix, 10 ms bins

unit_count = ev_unit_mx.shape[0]

# ----- for SUSTAINED - z-score, smooth and take 1st PC

# 1. sum all spikes for each pulse BEFORE z-scoring and PCA. 
# Works, but not in line with the previous analysis - 1 second smoothing is best

su_unit_mx_events = np.zeros([unit_count, len(sound_events)])
for i in range(unit_count):
    for j in range(len(sound_events)):
        su_unit_mx_events[i][j] = su_unit_mx[i][j*su_bin_count:(j+1)*su_bin_count].mean()
    
# smoothing
for j in range(unit_count):
    su_unit_mx_events[j] = smooth_rectangular(su_unit_mx_events[j], smooth_su_size)

# normalize
su_unit_mx_events_z = np.zeros([unit_count, len(sound_events)])
for j in range(unit_count):
    #su_unit_mx_events_z[j] = stats.zscore(su_unit_mx_events[j])
    su_unit_mx_events_z[j] = minmax_scale(su_unit_mx_events[j], feature_range=(0, 1), axis=0, copy=True)

# take first PC
su_pca = decomposition.PCA(n_components=2)
su_X   = su_pca.fit_transform(su_unit_mx_events_z.T)
su_resp = su_X[:, 0]  # PC1 scores. bin_size resolution

# OR just sum - not that cool
#su_resp = su_unit_mx_events_z.mean(axis=0)


# ----- for EVOKED - z-score, template match?, subtract sustained and sum / PCA

event_ids = {1: 'BGR', 2: 'TGT', 0: 'SIL', -1: 'NOI'}

# response profile matrix
psth_bins, psths_all = get_psth_matrix(snakemake.input[2], electrodes)
conditions = list(psths_all.keys())

# taking only the evoked profile part (important - this is not periodic!)
idx_s = int(psth_bins.shape[0]/2)
idx_e = idx_s + ev_bin_count # int(np.ceil(idx_s/2))

ev_unit_mx_events = np.zeros([len(sound_events), unit_count])
for i in range(unit_count):

    # subtract instantaneous FR of sustained part of that unit - not nice
    #idxs_del = ((np.arange(len(sound_events))+1)*su_bin_count)-1
    #su_resp_inst_unit = np.delete(su_unit_mx[i], idxs_del)
    #ev_unit_mx[i] = su_resp_inst_unit[:len(ev_unit_mx[i])]

    ev_unit_mx[i] = stats.zscore(ev_unit_mx[i])

    # subtract instantaneous sustained response multiplied by eigenvalue of that unit - not nice
    #su_inst = su_pca.components_[0][i] * su_resp  # subtract sustained part per unit
    #su_inst = np.repeat(su_inst.T, ev_bin_count)[:len(ev_unit_mx[i])]
    #ev_unit_mx[i] = ev_unit_mx[i] - su_inst

    for j in range(len(sound_events)):  # for each pulse do template matching via dot product
        idx_ev_mx = j*ev_bin_count
        #su_inst_unit = su_pca.components_[0][i] * su_resp[j]
        resp = ev_unit_mx[i][idx_ev_mx:idx_ev_mx+ev_bin_count] #- su_inst_unit
        if len(resp) == ev_bin_count:
            cond = event_ids[int(sound_events[j][1])]
            resp = np.dot(resp, psths_all[cond][:, idx_s:idx_e][i])  # evoked part only!
            ev_unit_mx_events[j][i] = resp

    #ev_unit_mx_events[:, i] = smooth_gaussian(ev_unit_mx_events[:, i], 25)

# take first PC - here its not cool
# ev_pca = decomposition.PCA(n_components=2)
# ev_X   = ev_pca.fit_transform(ev_unit_mx_events)
# ev_resp = ev_X[:, 0]  # PC1 scores. bin_size resolution

# OR just sum
ev_resp = ev_unit_mx_events.mean(axis=1)

# smooth evoked?
ev_resp = smooth_rectangular(ev_resp, 12)

# ----- finally dump everything ----------------

resp_manifold = np.column_stack([ev_resp, su_resp])

with h5py.File(snakemake.output[0], 'w') as f:
    f.create_dataset('unit_mx_proc_ev', data=ev_unit_mx_events)  # transformed activity matrix, event sampling
    f.create_dataset('unit_mx_proc_su', data=su_unit_mx_events.T)  # transformed activity matrix, event sampling
    f.create_dataset('evoked_response', data=ev_resp)
    f.create_dataset('sustained_response', data=su_resp)
    f.create_dataset('response_manifold', data=resp_manifold)

    for k, cond in enumerate(conditions):
        f.create_group(cond)  # sound conditions - BGR, TGT etc.
        f[cond].create_dataset('psth_mx', data=psths_all[cond])
        # f[cond].create_dataset('CCR_mx', data=CCR_mx_all[cond])