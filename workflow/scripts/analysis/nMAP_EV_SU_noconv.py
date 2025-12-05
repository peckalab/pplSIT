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

# 1. -------------- loading datasets -----------------------

# reading data
with h5py.File(snakemake.input[0], 'r') as f:
    sound_events = np.array(f['processed']['sound_events'])
    tl = np.array(f['processed']['timeline'])
    tgt_mx = np.array(f['processed']['target_matrix'])
    cfg = json.loads(f['processed'].attrs['parameters'])

# building unit activity matrix - shape (n_units, n_bins), 10 ms bins
bins, spikes_uxt = unit_activity_matrix(snakemake.input[0], snakemake.input[1], electrodes)

# reading PSTHs
event_ids = {1: 'BGR', 2: 'TGT', 0: 'SIL', -1: 'NOI'}
psth_bins, psths_all = get_psth_matrix(snakemake.input[2], electrodes)


# 2. -------------- computing population responses -----------------------

stim_onsets_bins = np.arange(len(sound_events))*25  # shape (n_trials,), integer bin index of each stimulus onset
stim_labels = sound_events[:, 1]  # per-trial labels (Target/Background/NoStim)
bin_size_s = 0.01
evk_ms = [0, 130]    # Evoked: 10–30 ms
sus_ms = [120, 250]  # Sustained: 120–250 ms

U, T = spikes_uxt.shape

# Convert ms windows to bin offsets
def ms2bins(ms): return int(round(ms / (bin_size_s * 1000.0)))
evk_b   = ( ms2bins(evk_ms[0]), ms2bins(evk_ms[1]) )
sus_b   = ( ms2bins(sus_ms[0]), ms2bins(sus_ms[1]) )

n_trials = len(stim_onsets_bins)

# Helper to extract counts in a relative bin window
def window_counts(rel_window, do_template_matching=True):
    start_off, end_off = rel_window
    # inclusive-exclusive indexing
    counts = np.zeros((U, n_trials), dtype=float)
    for ti, t0 in enumerate(stim_onsets_bins):
        s = t0 + start_off
        e = t0 + end_off
        if s < 0 or e > T:  # guard edges
            counts[:, ti] = np.nan
        else:
            if do_template_matching:  # align by PSTH template
                cond = event_ids[int(sound_events[ti][1])]
                unit_resp = spikes_uxt[:, s:e]
                mid_bin = psths_all[cond].shape[1]//2
                psth_template = psths_all[cond][:, mid_bin + start_off:mid_bin + end_off]
                counts[:, ti] = np.sum(unit_resp * psth_template, axis=1) / psth_template.sum(axis=1)
            else:  # sum raw counts across bins in window
                counts[:, ti] = spikes_uxt[:, s:e].sum(axis=1)
            
    return counts  # (n_units, n_trials)

# Raw counts per unit × trial for each window
counts_evk = window_counts(evk_b, do_template_matching=True)   # evoked
counts_sus = window_counts(sus_b, do_template_matching=False)   # sustained

# Variance-stabilize (Anscombe) to reduce Poisson sparsity issues
def anscombe(x): return 2.0 * np.sqrt(x + 3.0/8.0)
evk_a = anscombe(counts_evk)
sus_a = anscombe(counts_sus)

# Build each unit's baseline across all trials (drop NaNs)
# You can switch this to use only No-stimulus trials if you prefer.
#mu_pre = np.nanmean(pre_a, axis=1, keepdims=True)         # (n_units, 1)
#sd_pre = np.nanstd (pre_a, axis=1, ddof=1, keepdims=True) # (n_units, 1)
#sd_pre = np.where(sd_pre < 1e-6, 1e-6, sd_pre)            # avoid divide-by-zero

# Per-unit z-scores for evoked/sustained
#evoked_z    = (evk_a - mu_pre) / sd_pre                   # (n_units, n_trials)
#sustained_z = (sus_a - mu_pre) / sd_pre

# or z-score relative to the states
evoked_z    = (evk_a - np.nanmean(evk_a, axis=1, keepdims=True)) / np.nanstd(evk_a, axis=1, keepdims=True)
sustained_z = (sus_a - np.nanmean(sus_a, axis=1, keepdims=True)) / np.nanstd(sus_a, axis=1, keepdims=True)

# Population = mean across units (ignoring NaNs)
#evoked_pop_z    = np.nanmean(evoked_z, axis=0)            # (n_trials,)
#sustained_pop_z = np.nanmean(sustained_z, axis=0)

# or PCA
evoked_z    = np.nan_to_num(evoked_z, nan=0.0)
sustained_z = np.nan_to_num(sustained_z, nan=0.0)
EV_PC1 = decomposition.PCA(n_components=2).fit_transform(evoked_z.T)[:, 0]
SU_PC1 = decomposition.PCA(n_components=2).fit_transform(sustained_z.T)[:, 0]

# smoothing - large windows, but...
EV_PC1 = smooth_rectangular(EV_PC1, smooth_ev_size)
SU_PC1 = smooth_rectangular(SU_PC1, smooth_su_size)


# fix sign for evoked, also sustained - based only on the fact that no stim periods
# have usually more running. So SU values should be higher
idxs_bgr_ev = np.where(stim_labels == 1)[0]
idxs_sil_ev = np.where(stim_labels == 0)[0]
if np.nanmean(EV_PC1[idxs_bgr_ev]) < np.nanmean(EV_PC1[idxs_sil_ev]):
    EV_PC1 *= -1
if np.nanmean(SU_PC1[idxs_bgr_ev]) > np.nanmean(SU_PC1[idxs_sil_ev]):
    SU_PC1 *= -1
    
# 2) Center on no-stim / background
mu_ns = np.nanmean(EV_PC1[idxs_sil_ev])
EV_PC1_c = EV_PC1 - mu_ns
mu_bg = np.nanmean(SU_PC1[idxs_bgr_ev])
SU_PC1_c = SU_PC1 - mu_bg

ev_resp = EV_PC1_c
su_resp = SU_PC1_c


# 3. -------------- response manifold -----------------------

resp_manifold = np.column_stack([ev_resp, su_resp])

# finally dump everything
with h5py.File(snakemake.output[0], 'w') as f:
    f.create_dataset('unit_mx_proc_ev', data=evoked_z.T)  # time x units, transformed activity matrix, event sampling
    f.create_dataset('unit_mx_proc_su', data=sustained_z.T)  # time x units, transformed activity matrix, event sampling
    f.create_dataset('evoked_response', data=ev_resp)
    f.create_dataset('sustained_response', data=su_resp)
    f.create_dataset('response_manifold', data=resp_manifold)






# -------------------------------------
# ----------- OLD WAY -----------------
# -------------------------------------

# with h5py.File(snakemake.input[1], 'r') as f:
#     ev_bin_count = int(np.array(f['ev_bin_count']))
#     ev_periods = np.array(f['ev_periods'])
#     ev_bins = np.array(f['ev_bins'])
#     ev_unit_mx = np.array(f['ev_unit_mx'])  # original activity matrix, 10 ms bins

#     su_bin_count = int(np.array(f['su_bin_count']))
#     su_periods = np.array(f['su_periods'])
#     su_bins = np.array(f['su_bins'])
#     su_unit_mx = np.array(f['su_unit_mx'])  # original activity matrix, 10 ms bins

# unit_count = ev_unit_mx.shape[0]


# # ----- for SUSTAINED - z-score, smooth and take 1st PC

# # 1. sum all spikes for each pulse BEFORE z-scoring and PCA. 
# # Works, but not in line with the previous analysis - 1 second smoothing is best

# su_unit_mx_events = np.zeros([unit_count, len(sound_events)])
# for i in range(unit_count):
#     for j in range(len(sound_events)):
#         su_unit_mx_events[i][j] = su_unit_mx[i][j*su_bin_count:(j+1)*su_bin_count].mean()
    
# # smoothing
# for j in range(unit_count):
#     su_unit_mx_events[j] = smooth_rectangular(su_unit_mx_events[j], smooth_su_size)

# # normalize
# su_unit_mx_events_z = np.zeros([unit_count, len(sound_events)])
# for j in range(unit_count):
#     #su_unit_mx_events_z[j] = stats.zscore(su_unit_mx_events[j])
#     su_unit_mx_events_z[j] = minmax_scale(su_unit_mx_events[j], feature_range=(0, 1), axis=0, copy=True)

# # take first PC
# su_pca = decomposition.PCA(n_components=2)
# su_X   = su_pca.fit_transform(su_unit_mx_events_z.T)
# su_resp = su_X[:, 0]  # PC1 scores. bin_size resolution

# # 2. Sum spikes after PCA
# #for j in range(unit_count):
# #    su_unit_mx[j] = smooth_rectangular(su_unit_mx[j], snakemake.config['nMAP_EV_SU']['k_width'])
# #     su_unit_mx[j] = stats.zscore(su_unit_mx[j])

# # su_pca = decomposition.PCA(n_components=2)
# # su_X   = su_pca.fit_transform(su_unit_mx.T)
# # su_resp = su_X[:, 0]  # PC1 scores. bin_size resolution
# # su_resp = su_resp.reshape([int(su_resp.shape[0]/su_bin_count), su_bin_count]).mean(axis=1)  # sound events resolution

# # smooth again the final thing
# #su_resp = smooth_gaussian(su_resp, snakemake.config['nMAP_EV_SU']['smooth_k'])

# # TODO: resolve why smoothing AFTER z-scoring/PCA, not BEFORE.
# # try to smooth before, diff kernel sizes


# # ----- for EVOKED - z-score, using just the on-response window + PCA

# ev_unit_mx_events = np.zeros([unit_count, len(sound_events)])
# for i in range(unit_count):
#     for j in range(len(sound_events)):
#         unit_ev_resp = ev_unit_mx[i][j*ev_bin_count + 1:j*ev_bin_count + 3].mean()
#         if not np.isnan(unit_ev_resp):
#             ev_unit_mx_events[i][j] = unit_ev_resp
    
# # z-score / smoothing
# for j in range(unit_count):
#     ev_unit_mx_events[j] = smooth_rectangular(ev_unit_mx_events[j], smooth_ev_size)
#     ev_unit_mx_events[j] = ev_unit_mx_events[j] - su_unit_mx_events[j]  # subtract sustained
#     #ev_unit_mx_events[j] = stats.zscore(ev_unit_mx_events[j])
#     ev_unit_mx_events[j] = minmax_scale(ev_unit_mx_events[j], feature_range=(0, 1), axis=0, copy=True)

# # PCA:
# ev_pca = decomposition.PCA(n_components=2)
# ev_X   = ev_pca.fit_transform(ev_unit_mx_events.T)
# ev_resp = ev_X[:, 0]  # PC1 scores. sound events resolution

# # or just a sum - works better
# #ev_resp = ev_unit_mx_events.mean(axis=0)
# ev_resp = smooth_gaussian(ev_resp, smooth_ev_size)

# # TODO:
# # - try PCA on evoked
# # - try diff smoothing for sustained?
# # - try smoothing on evoked too

# # 3. -------------- response manifold -----------------------

# resp_manifold = np.column_stack([ev_resp, su_resp])

# # finally dump everything
# with h5py.File(snakemake.output[0], 'w') as f:
#     f.create_dataset('unit_mx_proc_ev', data=ev_unit_mx_events.T)  # transformed activity matrix, event sampling
#     f.create_dataset('unit_mx_proc_su', data=su_unit_mx_events.T)  # transformed activity matrix, event sampling
#     f.create_dataset('evoked_response', data=ev_resp)
#     f.create_dataset('sustained_response', data=su_resp)
#     f.create_dataset('response_manifold', data=resp_manifold)


# ------------------------------------------
# ----------- VERY OLD WAY -----------------
# ------------------------------------------

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