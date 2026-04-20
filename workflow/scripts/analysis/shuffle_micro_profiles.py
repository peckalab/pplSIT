import os
import sys

# prevent nested multithreading oversubscription
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import h5py
import numpy as np
from joblib import Parallel, delayed

# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.psth import staple_pulsetrain, staple_spike_times, get_spike_counts
from utils.events import get_event_periods
from utils.neurosuite import get_unit_names_sorted


def get_shuffled_with_rng(spiketrain, rng):
    """
    Shuffle spike times while preserving the ISI distribution.
    """
    if spiketrain.size <= 1:
        return spiketrain.copy()

    isis = np.diff(spiketrain).copy()
    rng.shuffle(isis)
    return np.concatenate(([spiketrain[0]], spiketrain[0] + np.cumsum(isis)))


def compute_shuffled_metrics_rng(strain, event_times, offset, bin_count, iter_count, rng):
    """
    Same output format as legacy compute_shuffled_metrics:
    rows = [bins[:-1], mean, std, p2.5, p97.5]
    """
    psth_shuffled = np.zeros((iter_count, bin_count - 1), dtype=np.float32)

    bins = None
    for i in range(iter_count):
        shuffled = get_shuffled_with_rng(strain, rng)
        bins, psth = get_spike_counts(shuffled, event_times, hw=offset, bin_count=bin_count)
        psth_shuffled[i] = psth.astype(np.float32)

    mean = psth_shuffled.mean(axis=0)
    std = psth_shuffled.std(axis=0)
    lo = np.percentile(psth_shuffled, 2.5, axis=0)
    hi = np.percentile(psth_shuffled, 97.5, axis=0)

    stats = np.vstack([
        bins[:-1],
        mean,
        std,
        lo,
        hi,
    ]).astype(np.float32)

    return stats


def compute_unit_event_shuffled(
    j, k, event_name,
    adjusted_pulses,
    periods,
    spike_times_unit,
    hw, bc, iter_count,
    seed_base=12345
):
    """
    Compute shuffled PSTH stats for one (event, unit) pair.
    """
    rng = np.random.default_rng(seed_base + j * 1_000_000 + k)

    # staple spike train within selected periods
    strain = staple_spike_times(spike_times_unit, periods, mode='sequence')
    strain = np.concatenate(strain) if len(strain) > 0 else np.array([], dtype=np.float64)

    # handle empty spike train safely
    if strain.size == 0 or adjusted_pulses.size == 0:
        bins = np.linspace(-hw, hw, bc)[:-1].astype(np.float32)
        zeros = np.zeros(bc - 1, dtype=np.float32)
        stats = np.vstack([bins, zeros, zeros, zeros, zeros]).astype(np.float32)
        return j, k, stats

    stats = compute_shuffled_metrics_rng(
        strain=strain,
        event_times=adjusted_pulses,
        offset=hw,
        bin_count=bc,
        iter_count=iter_count,
        rng=rng,
    )

    return j, k, stats


# ---------------------------
# configuration
# ---------------------------
hw = snakemake.config['psth']['bootstrap']['latency']
bc = snakemake.config['psth']['bootstrap']['bin_count']
iter_count = snakemake.config['psth']['bootstrap']['boot_iter_count']

event_types = [0, 1, 2, -1]  # SIL, BGR, TGT, NOI
event_type_names = ['SIL', 'BGR', 'TGT', 'NOI']

n_jobs = int(getattr(snakemake, "threads", 1))

# ---------------------------
# read input data
# ---------------------------
with h5py.File(snakemake.input[0], 'r') as f:
    tl = np.array(f['processed']['timeline'])
    sound_events = np.array(f['processed']['sound_events'])

spike_times = {}
with h5py.File(snakemake.input[1], 'r') as f:
    unit_names = get_unit_names_sorted([name for name in f])
    for unit_name in unit_names:
        spike_times[unit_name] = np.sort(np.array(f[unit_name]['spike_times']))

# ---------------------------
# precompute event-specific pulse trains and periods
# ---------------------------
event_info = []
for event_id, event_name in zip(event_types, event_type_names):
    pulses = sound_events[sound_events[:, 1] == event_id][:, 0]
    periods = get_event_periods(tl, event_id)

    adjusted_pulses = staple_pulsetrain(pulses, periods)

    event_info.append({
        "event_id": event_id,
        "event_name": event_name,
        "periods": periods,
        "adjusted_pulses": adjusted_pulses,
    })

# ---------------------------
# allocate output
# ---------------------------
J = len(event_info)
U = len(unit_names)

shuffled = np.zeros((J, U, 5, bc - 1), dtype=np.float32)

# ---------------------------
# build tasks
# ---------------------------
tasks = []
for j, ev in enumerate(event_info):
    for k, unit_name in enumerate(unit_names):
        tasks.append((
            j,
            k,
            ev["event_name"],
            ev["adjusted_pulses"],
            ev["periods"],
            spike_times[unit_name],
        ))

n_jobs = min(n_jobs, len(tasks)) if len(tasks) > 0 else 1

# ---------------------------
# parallel execution
# ---------------------------
results = Parallel(n_jobs=n_jobs, backend="loky", batch_size=1)(
    delayed(compute_unit_event_shuffled)(
        j=j,
        k=k,
        event_name=event_name,
        adjusted_pulses=adjusted_pulses,
        periods=periods,
        spike_times_unit=spike_times_unit,
        hw=hw,
        bc=bc,
        iter_count=iter_count,
    )
    for (j, k, event_name, adjusted_pulses, periods, spike_times_unit) in tasks
)

for j, k, stats in results:
    shuffled[j, k] = stats

# ---------------------------
# save to H5
# ---------------------------
with h5py.File(snakemake.output[0], 'w') as f:
    for i, event_name in enumerate(event_type_names):
        grp_ev = f.create_group(event_name)

        for j, unit_name in enumerate(unit_names):
            grp_unit = grp_ev.create_group(unit_name)
            grp_unit.create_dataset('shuffled', data=shuffled[i, j])


# import os, sys
# import h5py
# import numpy as np

# # import util functions from utils module
# parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
# sys.path.append(os.getcwd())
# sys.path.append(parent_dir)

# from utils.psth import compute_shuffled_metrics, staple_pulsetrain, staple_spike_times
# from utils.events import get_event_periods
# from utils.neurosuite import get_unit_names_sorted

# # some configs
# hw = snakemake.config['psth']['bootstrap']['latency']
# bc = snakemake.config['psth']['bootstrap']['bin_count']
# iter_count = snakemake.config['psth']['bootstrap']['boot_iter_count']
# event_types = [0, 1, 2, -1]  # SIL, BGR, TGT, NOI - order matters
# event_type_names = ['SIL', 'BGR', 'TGT', 'NOI']  # SIL, BGR, TGT, NOI - order matters

# # reading timeline and spike times
# with h5py.File(snakemake.input[0], 'r') as f:
#     tl = np.array(f['processed']['timeline'])
#     sound_events = np.array(f['processed']['sound_events'])
    
# spike_times = {}
# with h5py.File(snakemake.input[1], 'r') as f:
#     unit_names = get_unit_names_sorted([name for name in f])
#     for unit_name in f:
#         spike_times[unit_name] = np.array(f[unit_name]['spike_times'])

# # compute shuffled PSTHs
# shuffled = np.zeros([len(event_types), len(unit_names), 5, bc-1])
# for j, event_id in enumerate(event_types):  # sound events
#     pulses  = sound_events[sound_events[:, 1] == event_id][:, 0]
#     periods = get_event_periods(tl, event_id)

#     # shrink all selected pulses as if there is no other periods.
#     # allows to do shuffling without effects from other periods
#     # where mean firing rates can be different
#     adjusted_pulses = staple_pulsetrain(pulses, periods)
    
#     for k, unit_name in enumerate(unit_names):
#         s_times = spike_times[unit_name]
#         strain = staple_spike_times(s_times, periods, mode='sequence')  # result is in periods!
#         strain = np.array([item for sublist in strain for item in sublist])  # flatten to one array
#         shuffled[j][k] = compute_shuffled_metrics(strain, adjusted_pulses, offset=hw, bin_count=bc, iter_count=iter_count)

# # save to H5
# with h5py.File(snakemake.output[0], 'w') as f:
#     for i, event_name in enumerate(event_type_names):
#         grp_ev = f.create_group(event_name)

#         for j, unit_name in enumerate(unit_names):
#             grp_unit = grp_ev.create_group(unit_name)
#             grp_unit.create_dataset('shuffled', data=shuffled[i][j])