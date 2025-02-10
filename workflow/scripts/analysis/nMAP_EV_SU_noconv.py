import h5py, os, sys, json
import numpy as np
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

# 2. Sum spikes after PCA
#for j in range(unit_count):
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


# ----- for EVOKED - z-score, using just the on-response window + PCA

ev_unit_mx_events = np.zeros([unit_count, len(sound_events)])
for i in range(unit_count):
    for j in range(len(sound_events)):
        unit_ev_resp = ev_unit_mx[i][j*ev_bin_count + 1:j*ev_bin_count + 3].mean()
        if not np.isnan(unit_ev_resp):
            ev_unit_mx_events[i][j] = unit_ev_resp
    
# z-score / smoothing
for j in range(unit_count):
    ev_unit_mx_events[j] = smooth_rectangular(ev_unit_mx_events[j], smooth_ev_size)
    ev_unit_mx_events[j] = ev_unit_mx_events[j] - su_unit_mx_events[j]  # subtract sustained
    #ev_unit_mx_events[j] = stats.zscore(ev_unit_mx_events[j])
    ev_unit_mx_events[j] = minmax_scale(ev_unit_mx_events[j], feature_range=(0, 1), axis=0, copy=True)

# PCA:
ev_pca = decomposition.PCA(n_components=2)
ev_X   = ev_pca.fit_transform(ev_unit_mx_events.T)
ev_resp = ev_X[:, 0]  # PC1 scores. sound events resolution

# or just a sum - works better
#ev_resp = ev_unit_mx_events.mean(axis=0)
ev_resp = smooth_gaussian(ev_resp, smooth_ev_size)

# TODO:
# - try PCA on evoked
# - try diff smoothing for sustained?
# - try smoothing on evoked too



# 3. -------------- response manifold -----------------------

resp_manifold = np.column_stack([ev_resp, su_resp])

# finally dump everything
with h5py.File(snakemake.output[0], 'w') as f:
    f.create_dataset('unit_mx_proc_ev', data=ev_unit_mx_events.T)  # transformed activity matrix, event sampling
    f.create_dataset('unit_mx_proc_su', data=su_unit_mx_events.T)  # transformed activity matrix, event sampling
    f.create_dataset('evoked_response', data=ev_resp)
    f.create_dataset('sustained_response', data=su_resp)
    f.create_dataset('response_manifold', data=resp_manifold)



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