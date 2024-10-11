import h5py, os, sys, json
import numpy as np
from scipy import stats, signal
from sklearn import decomposition


# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.population import unit_activity_matrix
from utils.psth import get_psth_matrix

def make_smooth(data, k_width):
    kernel  = signal.gaussian(k_width, std=(k_width) / 7.2)
    return np.convolve(data, kernel, 'same') / kernel.sum()


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

# 1. -------------- Compute EVOKED response -----------------------

# activity matrix
bin_size = 0.01
latency  = cfg['sound']['latency']  # seconds
bins_per_event = int(latency / bin_size)
bins_ev, unit_mx_ev = unit_activity_matrix(snakemake.input[0], snakemake.input[1], electrodes, bin_size=bin_size)
    
# z-score
for i in range(len(unit_mx_ev)):
    unit_mx_ev[i] = stats.zscore(unit_mx_ev[i])

# response profile matrix
psth_bins, psths_all = get_psth_matrix(snakemake.input[2], electrodes)
conditions = list(psths_all.keys())

# taking only the evoked profile part (important - this is not periodic!)
idx_s = int(psth_bins.shape[0]/2)
idx_e = idx_s + int(np.ceil(idx_s/2))

# compute CCR matrix
CCR_mx_all = {}

for k, cond in enumerate(conditions):
    CCR_mx = np.zeros(unit_mx_ev.shape)
    for unit_idx in range(unit_mx_ev.shape[0]):
        prof = psths_all[cond][:, idx_s:idx_e][unit_idx]  # evoked part only!
        #prof = psths_all[cond][unit_idx]  # all
        spks = unit_mx_ev[unit_idx]
        CCR_mx[unit_idx] = signal.correlate(spks, prof, mode='same')
        
    CCR_mx_all[cond] = CCR_mx

# compute EVOKED response
evoked_resp = np.zeros(len(sound_events))
event_ids = {'BGR': 1, 'TGT': 2, 'SIL': 0, 'NOI': -1}

idx_peak = 6  # this number is very important - which phase of CCR to take
for k, cond in enumerate(conditions):
    sig = CCR_mx_all[cond].mean(axis=0)
    sig_cond = sig[idx_peak::bins_per_event]  

    ev_id = event_ids[cond]
    ev_idxs = np.where(sound_events[:, 1] == ev_id)[0]
    for idx in ev_idxs:
        evoked_resp[idx] = sig_cond[idx]

# smooth
evoked_resp_sm = make_smooth(evoked_resp, snakemake.config['nMAP_EV_SU']['k_width'])

# 2. -------------- Compute SUSTAINED response -----------------------

# activity matrix
bin_size = 0.125
bins_su, unit_mx_su = unit_activity_matrix(snakemake.input[0], snakemake.input[1], electrodes, bin_size=bin_size)

# z-score
for i in range(len(unit_mx_su)):
    unit_mx_su[i] = stats.zscore(unit_mx_su[i])

# new way - PCA on all units
su_mx = unit_mx_su[:, 1::2].T
pca = decomposition.PCA(n_components=2)
pca.fit(su_mx)
X = pca.transform(su_mx)
sustained_resp = X[:, 0]  # PC1 score

# smooth
sustained_resp_sm = make_smooth(sustained_resp, snakemake.config['nMAP_EV_SU']['k_width'])

# 3. -------------- response manifold -----------------------
resp_manifold = np.column_stack([evoked_resp_sm, sustained_resp_sm])

# finally dump everything
with h5py.File(snakemake.output[0], 'w') as f:
    f.create_dataset('unit_mx_ev', data=unit_mx_ev)
    f.create_dataset('unit_mx_su', data=unit_mx_su)
    f.create_dataset('bins_ev', data=bins_ev)
    f.create_dataset('bins_su', data=bins_su)
    f.create_dataset('evoked_resp', data=evoked_resp)
    f.create_dataset('evoked_resp_sm', data=evoked_resp_sm)
    f.create_dataset('sustained_resp', data=sustained_resp)
    f.create_dataset('sustained_resp_sm', data=sustained_resp_sm)
    f.create_dataset('response_manifold', data=resp_manifold)

    for k, cond in enumerate(conditions):
        f.create_group(cond)  # sound conditions - BGR, TGT etc.
        f[cond].create_dataset('psth_mx', data=psths_all[cond])
        f[cond].create_dataset('CCR_mx', data=CCR_mx_all[cond])
