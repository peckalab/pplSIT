import h5py, os, sys, json
import numpy as np
from scipy import stats, signal
from sklearn import decomposition


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

# reading data
with h5py.File(snakemake.input[0], 'r') as f:
    sound_events = np.array(f['processed']['sound_events'])
    tl = np.array(f['processed']['timeline'])
    tgt_mx = np.array(f['processed']['target_matrix'])
    cfg = json.loads(f['processed'].attrs['parameters'])

spike_times = {}
with h5py.File(snakemake.input[1], 'r') as f:
    unit_names = [x for x in f if int(x.split('-')[0]) in electrodes]
with h5py.File(snakemake.input[1], 'r') as f:
    for unit_name in unit_names:
        spike_times[unit_name]  = np.array(f[unit_name]['spike_times'])

# -------------- Compute population responses -----------------------

bin_size = snakemake.config['nMAP_EV_SU']['bin_size']
ev_su_lag = snakemake.config['nMAP_EV_SU']['ev_su_lag']
smooth_su_size = snakemake.config['nMAP_EV_SU']['smooth_su_size']
smooth_ev_size = snakemake.config['nMAP_EV_SU']['smooth_ev_size']

ev_bin_count = int(ev_su_lag/bin_size)
ev_times_beg = sound_events[:, 0]
ev_times_end = ev_times_beg + ev_su_lag
ev_periods = np.vstack([ev_times_beg, ev_times_end]).T

su_bin_count = int((cfg['sound']['latency'] - ev_su_lag)/bin_size)
su_times_beg = sound_events[:, 0] + ev_su_lag
su_times_end = su_times_beg + cfg['sound']['latency'] - ev_su_lag
su_periods = np.vstack([su_times_beg, su_times_end]).T

ev_bins = np.arange(0, np.diff(ev_periods, axis=1).sum(), bin_size)  # ignore last uneven bin
ev_unit_mx = np.zeros([len(unit_names), len(ev_bins)-1])
su_bins = np.arange(0, np.diff(su_periods, axis=1).sum(), bin_size)  # ignore last uneven bin
su_unit_mx = np.zeros([len(unit_names), len(su_bins)-1])

for k, unit_name in enumerate(unit_names):
    s_times = spike_times[unit_name]

    # shrink all spikes as if there is no other periods.
    ev_strain = staple_spike_times(s_times, ev_periods, mode='sequence')  # result is in periods!
    ev_strain = np.array([item for sublist in ev_strain for item in sublist])  # flatten to one array
    ev_unit_mx[k] = np.histogram(ev_strain, bins=ev_bins)[0]

    su_strain = staple_spike_times(s_times, su_periods, mode='sequence')  # result is in periods!
    su_strain = np.array([item for sublist in su_strain for item in sublist])  # flatten to one array
    su_unit_mx[k] = np.histogram(su_strain, bins=su_bins)[0]


# ----- for SUSTAINED - z-score, smooth and take 1st PC

# 1. sum all spikes for each pulse BEFORE z-scoring and PCA. 
# Works, but not in line with the previous analysis - 1 second smoothing is best

su_unit_mx_events = np.zeros([len(unit_names), len(sound_events)])
for i in range(len(unit_names)):
    for j in range(len(sound_events)):
        su_unit_mx_events[i][j] = su_unit_mx[i][j*su_bin_count:(j+1)*su_bin_count].mean()
    
# smoothing
for j in range(len(unit_names)):
    su_unit_mx_events[j] = smooth_rectangular(su_unit_mx_events[j], smooth_su_size)

# z-score
su_unit_mx_events_z = np.zeros([len(unit_names), len(sound_events)])
for j in range(len(unit_names)):
    su_unit_mx_events_z[j] = stats.zscore(su_unit_mx_events[j])

# take first PC
su_pca = decomposition.PCA(n_components=2)
su_X   = su_pca.fit_transform(su_unit_mx_events_z.T)
su_resp = su_X[:, 0]  # PC1 scores. bin_size resolution

# 2. Sum spikes after PCA
#for j in range(len(unit_names)):
#    su_unit_mx[j] = smooth_rectangular(su_unit_mx[j], snakemake.config['nMAP_EV_SU']['k_width'])
#     su_unit_mx[j] = stats.zscore(su_unit_mx[j])

# su_pca = decomposition.PCA(n_components=2)
# su_X   = su_pca.fit_transform(su_unit_mx.T)
# su_resp = su_X[:, 0]  # PC1 scores. bin_size resolution
# su_resp = su_resp.reshape([int(su_resp.shape[0]/su_bin_count), su_bin_count]).mean(axis=1)  # sound events resolution

# smooth again the final thing
#su_resp = smooth_gaussian(su_resp, snakemake.config['nMAP_EV_SU']['smooth_k'])

# TODO: resolve why smoothing AFTER z-scoring/PCA, not BEFORE.
# try to smooth before, diff kernel sizes



# ----- for EVOKED - z-score, template match?, subtract sustained and sum / PCA


# OPTION 1. Using template matching
# event_ids = {1: 'BGR', 2: 'TGT', 0: 'SIL', -1: 'NOI'}

# # response profile matrix
# psth_bins, psths_all = get_psth_matrix(snakemake.input[2], electrodes)
# conditions = list(psths_all.keys())

# # taking only the evoked profile part (important - this is not periodic!)
# idx_s = int(psth_bins.shape[0]/2)
# idx_e = idx_s + ev_bin_count # int(np.ceil(idx_s/2))

# ev_resp_mx = np.zeros([len(sound_events), len(unit_names)])
# for i in range(len(unit_names)):

#     # subtract instantaneous FR of sustained part of that unit - not nice
#     #idxs_del = ((np.arange(len(sound_events))+1)*su_bin_count)-1
#     #su_resp_inst_unit = np.delete(su_unit_mx[i], idxs_del)
#     #ev_unit_mx[i] = su_resp_inst_unit[:len(ev_unit_mx[i])]

#     ev_unit_mx[i] = stats.zscore(ev_unit_mx[i])

#     # subtract instantaneous sustained response multiplied by eigenvalue of that unit - not nice
#     #su_inst = su_pca.components_[0][i] * su_resp  # subtract sustained part per unit
#     #su_inst = np.repeat(su_inst.T, ev_bin_count)[:len(ev_unit_mx[i])]
#     #ev_unit_mx[i] = ev_unit_mx[i] - su_inst

#     for j in range(len(sound_events)):  # for each pulse do template matching via dot product
#         idx_ev_mx = j*ev_bin_count
#         #su_inst_unit = su_pca.components_[0][i] * su_resp[j]
#         resp = ev_unit_mx[i][idx_ev_mx:idx_ev_mx+ev_bin_count] #- su_inst_unit
#         if len(resp) == ev_bin_count:
#             cond = event_ids[int(sound_events[j][1])]
#             resp = np.dot(resp, psths_all[cond][:, idx_s:idx_e][i])  # evoked part only!
#             ev_resp_mx[j][i] = resp
        
# OPTION 2. Using just the on-response window + PCA
ev_unit_mx_events = np.zeros([len(unit_names), len(sound_events)])
for i in range(len(unit_names)):
    for j in range(len(sound_events)):
        unit_ev_resp = ev_unit_mx[i][j*ev_bin_count + 1:j*ev_bin_count + 3].mean()
        if not np.isnan(unit_ev_resp):
            ev_unit_mx_events[i][j] = unit_ev_resp
    
# z-score / smoothing
for j in range(len(unit_names)):
    ev_unit_mx_events[j] = smooth_rectangular(ev_unit_mx_events[j], 4)
    ev_unit_mx_events[j] = ev_unit_mx_events[j] - su_unit_mx_events[j]  # subtract sustained
    ev_unit_mx_events[j] = stats.zscore(ev_unit_mx_events[j])

# PCA:
ev_pca = decomposition.PCA(n_components=2)
ev_X   = ev_pca.fit_transform(ev_unit_mx_events.T)
ev_resp = ev_X[:, 0]  # PC1 scores. sound events resolution

# or just a sum - works better
# ev_resp = ev_resp_mx.sum(axis=1) / len(unit_names)
# ev_resp = smooth_gaussian(ev_resp, snakemake.config['nMAP_EV_SU']['smooth_k'])

# 3. -------------- response manifold -----------------------

resp_manifold = np.column_stack([ev_resp, su_resp])

# finally dump everything
with h5py.File(snakemake.output[0], 'w') as f:
    f.create_dataset('unit_mx_ev', data=ev_unit_mx)  # original activity matrix, 10 ms bins
    f.create_dataset('unit_mx_su', data=su_unit_mx)  # original activity matrix, 10 ms bins
    f.create_dataset('unit_mx_proc_ev', data=ev_unit_mx_events.T)  # transformed activity matrix, event sampling
    f.create_dataset('unit_mx_proc_su', data=su_unit_mx_events.T)  # transformed activity matrix, event sampling
    f.create_dataset('bins_ev', data=ev_bins)
    f.create_dataset('bins_su', data=su_bins)
    f.create_dataset('evoked_response', data=ev_resp)
    f.create_dataset('sustained_response', data=su_resp)
    f.create_dataset('response_manifold', data=resp_manifold)

    # for k, cond in enumerate(conditions):
    #     f.create_group(cond)  # sound conditions - BGR, TGT etc.
    #     f[cond].create_dataset('psth_mx', data=psths_all[cond])
    #     f[cond].create_dataset('CCR_mx', data=CCR_mx_all[cond])




# ----------- OLD WAY -----------------
# # 1. -------------- Compute EVOKED response -----------------------

# # activity matrix
# bin_size = 0.01
# latency  = cfg['sound']['latency']  # seconds
# bins_per_event = int(latency / bin_size)
# bins_ev, unit_mx_ev = unit_activity_matrix(snakemake.input[0], snakemake.input[1], electrodes, bin_size=bin_size)
    
# # z-score
# for i in range(len(unit_mx_ev)):
#     unit_mx_ev[i] = stats.zscore(unit_mx_ev[i])

# # response profile matrix
# psth_bins, psths_all = get_psth_matrix(snakemake.input[2], electrodes)
# conditions = list(psths_all.keys())

# # taking only the evoked profile part (important - this is not periodic!)
# idx_s = int(psth_bins.shape[0]/2)
# idx_e = idx_s + int(np.ceil(idx_s/2))

# # compute CCR matrix
# CCR_mx_all = {}

# for k, cond in enumerate(conditions):
#     CCR_mx = np.zeros(unit_mx_ev.shape)
#     for unit_idx in range(unit_mx_ev.shape[0]):
#         prof = psths_all[cond][:, idx_s:idx_e][unit_idx]  # evoked part only!
#         #prof = psths_all[cond][unit_idx]  # all
#         spks = unit_mx_ev[unit_idx]
#         CCR_mx[unit_idx] = signal.correlate(spks, prof, mode='same')
        
#     CCR_mx_all[cond] = CCR_mx

# # compute EVOKED response
# evoked_resp = np.zeros(len(sound_events))
# event_ids = {'BGR': 1, 'TGT': 2, 'SIL': 0, 'NOI': -1}

# idx_peak = 6  # this number is very important - which phase of CCR to take
# for k, cond in enumerate(conditions):
#     sig = CCR_mx_all[cond].mean(axis=0)
#     sig_cond = sig[idx_peak::bins_per_event]  

#     ev_id = event_ids[cond]
#     ev_idxs = np.where(sound_events[:, 1] == ev_id)[0]
#     for idx in ev_idxs:
#         evoked_resp[idx] = sig_cond[idx]

# # smooth
# evoked_resp_sm = make_smooth(evoked_resp, snakemake.config['nMAP_EV_SU']['k_width'])

# # 2. -------------- Compute SUSTAINED response -----------------------

# # activity matrix
# bin_size = 0.125
# bins_su, unit_mx_su = unit_activity_matrix(snakemake.input[0], snakemake.input[1], electrodes, bin_size=bin_size)

# # z-score
# for i in range(len(unit_mx_su)):
#     unit_mx_su[i] = stats.zscore(unit_mx_su[i])

# # new way - PCA on all units
# su_mx = unit_mx_su[:, 1::2].T
# pca = decomposition.PCA(n_components=2)
# pca.fit(su_mx)
# X = pca.transform(su_mx)
# sustained_resp = X[:, 0]  # PC1 score

# # smooth
# sustained_resp_sm = make_smooth(sustained_resp, snakemake.config['nMAP_EV_SU']['k_width'])