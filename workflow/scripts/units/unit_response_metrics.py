import h5py, json, os, sys
import numpy as np

# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.spiketrain import unit_response_metrics_per_condition
from utils.neurosuite import get_unit_names_sorted


# reading data
with h5py.File(snakemake.input[0], 'r') as f:
    sound_events = np.array(f['processed']['sound_events'])
    tl = np.array(f['processed']['timeline'])
    cfg = json.loads(f['processed'].attrs['parameters'])

spike_times = {}
with h5py.File(snakemake.input[1], 'r') as f:
    units_to_plot = get_unit_names_sorted([name for name in f])
    for unit_name in units_to_plot:
        spike_times[unit_name] = np.array(f[unit_name]['spike_times'])

# for the moment make conditions manually - 
#state_names = snakemake.config['units']['response_metric_conditions']
response_metric_conditions = ['idxs_tgt_sta_succ', 'idxs_bgr_sta', 'idxs_bgr_run', 'idxs_sil_sta', 'idxs_sil_run']
idxs_states = {}
with h5py.File(snakemake.input[2], 'r') as f:
    for st_name in response_metric_conditions:
        idxs_states[st_name] = np.array(f[st_name])

# add generic "overlapping" conditions
idxs_bgr_ev = np.where(sound_events[:, 1] == 1)[0]
idxs_tgt_ev = np.where(sound_events[:, 1] == 2)[0]
idxs_sil_ev = np.where(sound_events[:, 1] == 0)[0]

idxs_states['idxs_stimulus_on'] = np.union1d(idxs_bgr_ev, idxs_tgt_ev)
idxs_states['idxs_stimulus_off'] = idxs_sil_ev

all_unit_metrics = {}
for unit_id, spiketrain in spike_times.items():
    all_unit_metrics[unit_id] = unit_response_metrics_per_condition(spike_times[unit_id], sound_events[:, 0], idxs_states)
    
# convert to condition-first dict + matrices
results = {}
for condition_id in idxs_states.keys():
    cond_mx = np.zeros([len(spike_times), 14])  # 14 metrics should be
    for i, (unit_id, metrics) in enumerate(all_unit_metrics.items()):
        if condition_id in metrics:
            cond_mx[i] = np.array(list(metrics[condition_id].values()))

    results[condition_id] = cond_mx.copy()

# save to H5
with h5py.File(snakemake.output[0], 'w') as f:
    for condition_id, metric_mx in results.items():
        grp = f.create_group(condition_id)
        ds = grp.create_dataset('metrics', data=metric_mx)
        ds.attrs['units_order'] = ','.join([key for key in spike_times])
        ds.attrs['metric_names'] = 'mean_rate, ptp, rms, rms_norm, latency_ms, reliability, jitter_ms,\
              fano_factor, normalized_ptp, ptp_evoked, ptp_baseline, ptp_ratio'