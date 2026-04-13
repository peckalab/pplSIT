import os
import sys
import json

# prevent nested BLAS/OpenMP oversubscription
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

from utils.neurosuite import get_unit_names_sorted
from utils.states import get_state_as_periods


def get_phases(pulse_times, spk_times, offset=0.25):
    phases = []
    for p_time in pulse_times:
        selected = spk_times[(spk_times > p_time) & (spk_times < p_time + offset)]
        phases += [2 * np.pi * x / offset for x in selected - p_time]
    return np.array(phases, dtype=np.float64)


def get_shuffled_with_rng(spiketrain, rng):
    """
    Shuffle spike times preserving ISIs, using a local RNG.
    """
    spiketrain = np.asarray(spiketrain, dtype=np.float64)
    if spiketrain.size <= 1:
        return spiketrain.copy()

    isis = np.diff(spiketrain).copy()
    rng.shuffle(isis)
    return np.concatenate([[spiketrain[0]], spiketrain[0] + np.cumsum(isis)])


def compute_condition_unit_metrics(
    cond_name,
    idxs_to_phase,
    unit_id,
    spk_times,
    pulse_times,
    n_shuffles,
    ipi_for_shuffle,
    psth_boot_cache,
    psth_shuf_cache,
    offset=0.25,
    seed_base=12345,
):
    rng = np.random.default_rng(
        seed_base
        + abs(hash(cond_name)) % 1_000_000
        + abs(hash(unit_id)) % 1_000_000
    )

    # real spike phases
    phases = get_phases(pulse_times[idxs_to_phase], spk_times, offset=offset)

    if len(phases) < 10:
        result = {
            "MRL_real": 0.0,
            "MRLs_shuffled": np.zeros(n_shuffles, dtype=np.float32),
            "p_value": 1.0,
        }
        return cond_name, unit_id, result

    # real MRL
    MRL_real = float(np.abs(np.mean(np.exp(1j * phases))))

    # staple spikes / pulses for shuffling
    shift = 0.0
    spikes_adjusted = []
    pulses_adjusted = []

    for i in idxs_to_phase:
        if i + 1 >= len(pulse_times):
            continue

        selected = spk_times[(spk_times > pulse_times[i]) & (spk_times < pulse_times[i + 1])]
        spikes_adjusted += [x + shift for x in (selected - pulse_times[i])]
        pulses_adjusted.append(shift)
        shift += ipi_for_shuffle

    spikes_adjusted = np.asarray(spikes_adjusted, dtype=np.float64)
    pulses_adjusted = np.asarray(pulses_adjusted, dtype=np.float64)

    # shuffle controls
    MRLs_shuffled = np.zeros(n_shuffles, dtype=np.float32)

    if len(spikes_adjusted) < 2 or len(pulses_adjusted) == 0:
        MRLs_shuffled[:] = 0.0
    else:
        for s in range(n_shuffles):
            strain_shuf = get_shuffled_with_rng(spikes_adjusted, rng)
            phases_shuf = get_phases(pulses_adjusted, strain_shuf, offset=offset)

            if len(phases_shuf) == 0:
                MRLs_shuffled[s] = 0.0
            else:
                MRLs_shuffled[s] = np.abs(np.mean(np.exp(1j * phases_shuf)))

    p_value = float(np.mean(MRLs_shuffled >= MRL_real))

    result = {
        "MRL_real": MRL_real,
        "MRLs_shuffled": MRLs_shuffled,
        "p_value": p_value,
    }

    # optional PSTH deviation score for BGR / TGT
    if cond_name in ("bgr", "tgt"):
        psth_boot = psth_boot_cache[cond_name][unit_id]
        psth_shuf = psth_shuf_cache[cond_name][unit_id]

        psth = psth_boot[1]
        CI_low = psth_shuf[3]
        CI_high = psth_shuf[4]
        shuf_mean = float(psth_shuf[1].mean())

        if shuf_mean != 0:
            idxs_above = np.where(psth - CI_high > 0)[0]
            idxs_below = np.where(psth - CI_low < 0)[0]

            dev_score = (
                (psth - CI_high)[idxs_above].sum()
                + np.abs((psth - CI_low)[idxs_below].sum())
            ) / shuf_mean
        else:
            dev_score = 0.0

        result["psth_dev_score"] = float(dev_score)

    return cond_name, unit_id, result


# -------------------------
# config / paths
# -------------------------
n_shuffles = snakemake.config["sound_phase_lock"]["n_shuffles"]
ipi_for_shuffle = snakemake.config["sound_phase_lock"]["ipi_for_shuffle"]
n_jobs = int(getattr(snakemake, "threads", 1))

s_path = os.path.dirname(snakemake.input[0])
session = os.path.basename(s_path)

# -------------------------
# reading events and spiking data
# -------------------------
with h5py.File(snakemake.input[0], "r") as f:
    tl = np.array(f["processed"]["timeline"])
    sound_events = np.array(f["processed"]["sound_events"])
    cfg = json.loads(f["processed"].attrs["parameters"])
    tgt_mx = np.array(f["processed"]["target_matrix"])

spike_times = {}
with h5py.File(snakemake.input[1], "r") as f:
    unit_names = get_unit_names_sorted([name for name in f])
    for unit_name in unit_names:
        spike_times[unit_name] = np.sort(np.array(f[unit_name]["spike_times"]))

# -------------------------
# build event indices
# -------------------------
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

# -------------------------
# build conditions
# -------------------------
bgr_sta_mx, idxs_bgr_sta_ev = get_state_as_periods(s_path, "BGR", "STA", None, 4)
bgr_run_mx, idxs_bgr_run_ev = get_state_as_periods(s_path, "BGR", "RUN", None, 2, strip_l=1)

cond_idxs = {
    "bgr": idxs_bgr_ev,
    "tgt": idxs_tgt_ev,
    "bgr_sta": idxs_bgr_sta_ev,
    "bgr_run": idxs_bgr_run_ev,
}

ensembles_f = os.path.join(s_path, "analysis", "ensembles.h5")
if os.path.exists(ensembles_f):
    bgr_sta_al_mx, idxs_bgr_sta_al_ev = get_state_as_periods(s_path, "BGR", "STA", "AL", 4)
    bgr_sta_ph_mx, idxs_bgr_sta_ph_ev = get_state_as_periods(s_path, "BGR", "STA", "PH", 2, strip_l=1)

    cond_idxs["bgr_sta_al"] = idxs_bgr_sta_al_ev
    cond_idxs["bgr_sta_ph"] = idxs_bgr_sta_ph_ev

# -------------------------
# preload PSTH stats once
# -------------------------
psth_boot_cache = {"bgr": {}, "tgt": {}}
psth_shuf_cache = {"bgr": {}, "tgt": {}}

with h5py.File(snakemake.input[2], "r") as f_boot, h5py.File(snakemake.input[3], "r") as f_shuf:
    for cond_name in ("bgr", "tgt"):
        grp_name = cond_name.upper()
        for unit_id in unit_names:
            psth_boot_cache[cond_name][unit_id] = np.array(f_boot[f"{grp_name}/{unit_id}/profile_stats"])
            psth_shuf_cache[cond_name][unit_id] = np.array(f_shuf[f"{grp_name}/{unit_id}/shuffled"])

# -------------------------
# parallel tasks
# -------------------------
tasks = []
for cond_name, idxs_to_phase in cond_idxs.items():
    idxs_to_phase = np.asarray(idxs_to_phase, dtype=np.int64)
    for unit_id in unit_names:
        tasks.append((cond_name, idxs_to_phase, unit_id))

n_jobs = min(n_jobs, len(tasks)) if len(tasks) > 0 else 1

results = Parallel(n_jobs=n_jobs, backend="loky", batch_size=1)(
    delayed(compute_condition_unit_metrics)(
        cond_name=cond_name,
        idxs_to_phase=idxs_to_phase,
        unit_id=unit_id,
        spk_times=spike_times[unit_id],
        pulse_times=pulse_times,
        n_shuffles=n_shuffles,
        ipi_for_shuffle=ipi_for_shuffle,
        psth_boot_cache=psth_boot_cache,
        psth_shuf_cache=psth_shuf_cache,
    )
    for cond_name, idxs_to_phase, unit_id in tasks
)

# -------------------------
# collect results
# -------------------------
unit_MRLs = {cond_name: {} for cond_name in cond_idxs.keys()}

for cond_name, unit_id, result in results:
    unit_MRLs[cond_name][unit_id] = result

# optional concise summary
for cond_name in cond_idxs.keys():
    print(f"{session}: phase-lock done for condition {cond_name} ({len(unit_MRLs[cond_name])} units)")

# -------------------------
# dump to H5
# -------------------------
with h5py.File(snakemake.output[0], "w") as f:
    for condition, data in unit_MRLs.items():
        grp_cond = f.create_group(condition)

        for unit_id, MRLs in data.items():
            grp_unit = grp_cond.create_group(unit_id)
            grp_unit.create_dataset("MRL_real", data=MRLs["MRL_real"])
            grp_unit.create_dataset("MRLs_shuffled", data=MRLs["MRLs_shuffled"])
            grp_unit.create_dataset("p_value", data=MRLs["p_value"])

            if "psth_dev_score" in MRLs:
                grp_unit.create_dataset("psth_dev_score", data=MRLs["psth_dev_score"])