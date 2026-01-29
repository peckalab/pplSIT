import os, sys

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import h5py
import json
import numpy as np
import matplotlib.pyplot as plt


# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.psth import get_spike_counts
from utils.neurosuite import get_unit_names_sorted

from joblib import Parallel, delayed

def compute_unit_state(j, u, idxs_pool, unit_name,
                       sound_events, spike_times_unit,
                       hw, bc, iter_count, bin_size, seed_base=12345):
    rng = np.random.default_rng(seed_base + j * 1_000_000 + u)

    prof = np.empty((iter_count, bc - 1), dtype=np.float32)

    for k in range(iter_count):
        idxs_rand = rng.choice(idxs_pool, size=len(idxs_pool), replace=True)
        times_rand = sound_events[idxs_rand, 0]

        # jitter (do not mutate original)
        strain = spike_times_unit + ((rng.random(spike_times_unit.shape[0]) - 0.5) * bin_size)

        bins, psth = get_spike_counts(strain, times_rand, hw=hw, bin_count=bc)
        prof[k] = psth

    med = np.median(prof, axis=0)
    std = prof.std(axis=0)
    lo  = np.percentile(prof,  2.5, axis=0)
    hi  = np.percentile(prof, 97.5, axis=0)

    stats = np.vstack([bins[:-1], med, std, lo, hi]).astype(np.float32)
    return j, u, prof, stats

# configuration
hw = snakemake.config['psth']['bootstrap']['latency']
bc = snakemake.config['psth']['bootstrap']['bin_count']
iter_count = snakemake.config['psth']['bootstrap']['boot_iter_count']
bin_size   = hw/((bc-1)/2)


# reading events and spiking data
with h5py.File(snakemake.input[0], 'r') as f:
    tl = np.array(f['processed']['timeline'])
    sound_events = np.array(f['processed']['sound_events'])
    cfg = json.loads(f['processed'].attrs['parameters'])

spike_times = {}
with h5py.File(snakemake.input[1], 'r') as f:
    unit_names = get_unit_names_sorted([name for name in f])
    for unit_name in f:
        spike_times[unit_name] = np.sort(
            np.array(f[unit_name]['spike_times'])
        )

#event_types = [0, 1, 2, -1]  # SIL, BGR, TGT, NOI - order matters
#event_type_names = ['SIL', 'BGR', 'TGT', 'NOI']  # SIL, BGR, TGT, NOI - order matters

# basic states
event_idxs_pool = {
    'SIL': np.where(sound_events[:, 1] == 0)[0],
    'BGR': np.where(sound_events[:, 1] == 1)[0],
    'TGT': np.where(sound_events[:, 1] == 2)[0],
    'NOI': np.where(sound_events[:, 1] ==-1)[0]
}
# distractor basic states
distr_count = int(cfg['experiment']['distractor_islands'])
if cfg['sound']['sounds']['distractor1']['enabled'] and distr_count > 0:
    event_idxs_pool['DI1'] = np.where(sound_events[:, 1] == 3)[0]
if cfg['sound']['sounds']['distractor2']['enabled'] and distr_count > 0:
    event_idxs_pool['DI2'] = np.where(sound_events[:, 1] == 4)[0]

# complex states
state_names  = ['idxs_tgt_sta_succ', 'idxs_bgr_sta', 'idxs_bgr_run', 'idxs_sil_sta', 'idxs_sil_run']
# distractor complex states
if cfg['sound']['sounds']['distractor1']['enabled'] and distr_count > 0:
    state_names.append('idxs_di1_sta')
if cfg['sound']['sounds']['distractor2']['enabled'] and distr_count > 0:
    state_names.append('idxs_di2_sta')
# if distractor fail
if cfg['experiment']['distractor_fail']:
    state_names.append('idxs_dis_fail')

with h5py.File(snakemake.input[2], 'r') as f:
    for st_name in state_names:
        event_idxs_pool[st_name] = np.array(f[st_name])

# now do bootstrapping, using parallel processing
n_jobs = int(getattr(snakemake, "threads", 1))



J = len(event_idxs_pool)
U = len(unit_names)

profiles = np.zeros((J, U, iter_count, bc - 1), dtype=np.float32)
profile_stats = np.zeros((J, U, 5, bc - 1), dtype=np.float32)

tasks = []
for j, (state_id, idxs_pool) in enumerate(event_idxs_pool.items()):
    for u, unit_name in enumerate(unit_names):
        tasks.append((j, u, idxs_pool, unit_name))

# clamp workers to number of tasks (optional but safe)
n_jobs = min(n_jobs, len(tasks))

results = Parallel(n_jobs=n_jobs, backend="loky", batch_size=1)(
    delayed(compute_unit_state)(
        j, u, idxs_pool, unit_name,
        sound_events, spike_times[unit_name],
        hw, bc, iter_count, bin_size
    )
    for (j, u, idxs_pool, unit_name) in tasks
)

for j, u, prof, stats in results:
    profiles[j, u] = prof
    profile_stats[j, u] = stats

# save to H5
with h5py.File(snakemake.output[0], 'w') as f:
    for i, event_name in enumerate(event_idxs_pool.keys()):
        grp_ev = f.create_group(event_name)

        for j, unit_name in enumerate(unit_names):
            grp_unit = grp_ev.create_group(unit_name)
            grp_unit.create_dataset('profiles', data=profiles[i][j])
            grp_unit.create_dataset('profile_stats', data=profile_stats[i][j])
