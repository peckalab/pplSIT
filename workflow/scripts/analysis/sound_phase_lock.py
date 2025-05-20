import os, sys
import h5py, json
import numpy as np

# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.neurosuite import get_unit_names_sorted
from utils.psth import get_shuffled
from utils.states import get_state_as_periods


def get_phases(pulse_times, spk_times, offset=0.25):
    phases = []
    for i, p_time in enumerate(pulse_times):
        selected = spk_times[(spk_times > p_time) & (spk_times < p_time + offset)]
        phases += [2 * np.pi * x/offset for x in selected - p_time]
    
    return np.array(phases)


# some configs
n_shuffles = snakemake.config['sound_phase_lock']['n_shuffles']
ipi_for_shuffle = snakemake.config['sound_phase_lock']['ipi_for_shuffle']
s_path  = os.path.dirname(snakemake.input[0])
session = os.path.basename(s_path)

# reading events and spiking data
with h5py.File(snakemake.input[0], 'r') as f:
    tl = np.array(f['processed']['timeline'])
    sound_events = np.array(f['processed']['sound_events'])
    cfg = json.loads(f['processed'].attrs['parameters'])
    tgt_mx = np.array(f['processed']['target_matrix'])

spike_times = {}
with h5py.File(snakemake.input[1], 'r') as f:
    unit_names = get_unit_names_sorted([name for name in f])
    for unit_name in f:
        spike_times[unit_name] = np.array(f[unit_name]['spike_times'])

# indices for diff conditions
x_pos_ev = tl[sound_events[:, 2].astype(np.int32)][:, 1]
y_pos_ev = tl[sound_events[:, 2].astype(np.int32)][:, 2]
speed_ev = tl[sound_events[:, 2].astype(np.int32)][:, 3]

pulse_times = sound_events[:, 0]

idxs_sta_ev = np.where(speed_ev < 0.04)[0]
idxs_run_ev = np.where(speed_ev > 0.04)[0]
idxs_bgr_ev = np.where(sound_events[:, 1] == 1)[0]
idxs_sil_ev = np.where(sound_events[:, 1] == 0)[0]
idxs_tgt_ev = np.where(sound_events[:, 1] == 2)[0]
idxs_di1_ev = np.where(sound_events[:, 1] == 3)[0]
idxs_di2_ev = np.where(sound_events[:, 1] == 4)[0]


# build conditions
bgr_sta_mx, idxs_bgr_sta_ev = get_state_as_periods(s_path, 'BGR', 'STA', None, 4)  # kind of correct way
bgr_run_mx, idxs_bgr_run_ev = get_state_as_periods(s_path, 'BGR', 'RUN', None, 2, strip_l=1)  # kind of correct way
cond_idxs  = {
    'bgr': idxs_bgr_ev,
    'tgt': idxs_tgt_ev,
    'bgr_sta': idxs_bgr_sta_ev,
    'bgr_run': idxs_bgr_run_ev
}

# if state ensembles exist, compute for them too
ensembles_f = os.path.join(s_path, 'analysis', 'ensembles.h5')
if os.path.exists(ensembles_f):
    bgr_sta_al_mx, idxs_bgr_sta_al_ev = get_state_as_periods(s_path, 'BGR', 'STA', 'AL', 4)
    bgr_sta_al_mx, idxs_bgr_sta_ph_ev = get_state_as_periods(s_path, 'BGR', 'STA', 'PH', 2, strip_l=1)

    cond_idxs['bgr_sta_al'] = idxs_bgr_sta_al_ev
    cond_idxs['bgr_sta_ph'] = idxs_bgr_sta_ph_ev

unit_MRLs = {}

# computing phases for diff conditions
for k, (cond_name, idxs_to_phase) in enumerate(cond_idxs.items()):
    MRLs_by_event_type = {}
    for j, (unit_id, spk_times) in enumerate(spike_times.items()):
        
        # get all spike phases
        phases = get_phases(pulse_times[idxs_to_phase], spk_times)
        
        # no spikes
        if len(phases) < 10:
            MRLs_by_event_type[unit_id] = { 
                "MRL_real": 0,
                "MRLs_shuffled": np.zeros(n_shuffles),
                "p_value": 1.0
            }
            continue
            
        # real MRL
        MRL_real = np.abs(np.mean(np.exp(1j * np.array(phases))))

        # staple spikes / pulses for shuffle controls
        shift = 0
        spikes_adjusted = []
        pulses_adjusted = []
        for i in idxs_to_phase:
            if i+1 >= len(pulse_times):
                continue

            selected = spk_times[(spk_times > pulse_times[i]) & (spk_times < pulse_times[i+1])]
            spikes_adjusted += [x + shift for x in selected - pulse_times[i]]
            pulses_adjusted += [shift]
            shift += ipi_for_shuffle

        # do shuffle controls
        MRLs_shuffled = []
        for _ in range(n_shuffles):
            strain_shuf = get_shuffled(spikes_adjusted)
            phases_shuf = get_phases(pulses_adjusted, strain_shuf)
            MRLs_shuffled.append(np.abs(np.mean(np.exp(1j * np.array(phases_shuf)))))
        MRLs_shuffled = np.array(MRLs_shuffled)

        # fraction of shuffled MRLs greater than real MRL
        p_value = np.mean(MRLs_shuffled >= MRL_real)

        MRLs_by_event_type[unit_id] = { 
            "MRL_real": MRL_real,
            "MRLs_shuffled": MRLs_shuffled,
            "p_value": p_value
        }
        print(f"{session}: unit {unit_id} phase lock done ({j} from {len(spike_times)}); ev type: {cond_name}")
        
    unit_MRLs[cond_name] = dict(MRLs_by_event_type)

# dump to H5
with h5py.File(snakemake.output[0], 'w') as f:
    for condition, data in unit_MRLs.items():
        grp_cond = f.create_group(condition)

        for unit_id, MRLs in data.items():
            grp_unit = grp_cond.create_group(unit_id)
            grp_unit.create_dataset('MRL_real', data=MRLs['MRL_real'])  # one number
            grp_unit.create_dataset('MRLs_shuffled', data=MRLs['MRLs_shuffled'])  # array of MRLs for each random shuffle
            grp_unit.create_dataset('p_value', data=MRLs['p_value'])  # one number
