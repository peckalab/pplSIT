import os, sys
import h5py
import numpy as np

# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.psth import compute_shuffled_metrics, staple_pulsetrain, staple_spike_times
from utils.events import get_event_periods


def compute_shuffle_psth_micro(meta_file, units_file, dst_unit_file, event_type, latency, bin_count):
    unit_name = os.path.basename(dst_unit_file).split('_')[0]
    with h5py.File(units_file, 'r') as f:
        s_times = np.array(f[unit_name]['spike_times'])

    with h5py.File(meta_file, 'r') as f:
        tl = np.array(f['processed']['timeline'])
        sound_events = np.array(f['processed']['sound_events'])

    pulses  = sound_events[sound_events[:, 1] == event_type][:, 0]
    periods = get_event_periods(tl, event_type)

    # first shrink all selected pulses as if there is no other periods
    adjusted_pulses = staple_pulsetrain(pulses, periods)

    # spikes during certain event type only
    strain = staple_spike_times(s_times, periods, mode='sequence')  # result is in periods!
    strain = np.array([item for sublist in strain for item in sublist])  # flatten to one array
    data = compute_shuffled_metrics(strain, adjusted_pulses, offset=latency, bin_count=bin_count)

    np.save(dst_unit_file, data)


# compute for silence
compute_shuffle_psth_micro(
    snakemake.input[0], 
    snakemake.input[1], 
    snakemake.output[0],
    0,  # silence event type
    snakemake.config['psth']['micro']['latency'],
    snakemake.config['psth']['micro']['bin_count']
)
