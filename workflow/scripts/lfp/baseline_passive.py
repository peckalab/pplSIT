import json
import os
import sys

import h5py
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.session import build_stimulus_catalog, detect_session_paradigm


def _decode_attr(value):
    if isinstance(value, bytes):
        return value.decode()
    return value


def _load_meta(meta_path):
    with h5py.File(meta_path, "r") as f:
        sound_events = np.array(f["processed"]["sound_events"])
        tl = np.array(f["processed"]["timeline"])
        cfg = json.loads(f["processed"].attrs["parameters"])
        paradigm = _decode_attr(f["processed"].attrs.get("session_paradigm", detect_session_paradigm(cfg)))
        stimulus_catalog_raw = _decode_attr(f["processed"].attrs.get("stimulus_catalog", ""))

    if stimulus_catalog_raw:
        stimulus_catalog = json.loads(stimulus_catalog_raw)
    else:
        stimulus_catalog = build_stimulus_catalog(cfg)

    return sound_events, tl, cfg, paradigm, stimulus_catalog


def _duration_sec(sound_id, stimulus_catalog):
    entry = stimulus_catalog.get(str(int(sound_id)), {})
    duration = entry.get("duration")
    if duration is None:
        return None
    duration = float(duration)
    if duration < 0:
        return None
    return duration


def _time_to_sample(time_s, fs_hz, n_samples):
    return int(np.clip(np.round(time_s * fs_hz), 0, n_samples))


def _mark_mask(mask, start_s, end_s, fs_hz):
    start_idx = _time_to_sample(start_s, fs_hz, mask.size)
    end_idx = _time_to_sample(end_s, fs_hz, mask.size)
    if end_idx > start_idx:
        mask[start_idx:end_idx] = True


def _resolve_event_end_time(sound_events, idx, stimulus_catalog, session_end_s):
    onset_s = float(sound_events[idx, 0])
    sound_id = int(sound_events[idx, 1])
    next_onset_s = session_end_s if idx + 1 >= len(sound_events) else float(sound_events[idx + 1, 0])
    next_onset_s = max(next_onset_s, onset_s)

    duration_s = _duration_sec(sound_id, stimulus_catalog)
    if duration_s is None:
        return next_onset_s

    return min(onset_s + duration_s, next_onset_s)


def _summarize_bootstrapped_means(boot_means, ci_percentile):
    boot_means = np.stack(boot_means, axis=0)
    lower_percentile = (100.0 - ci_percentile) / 2.0
    upper_percentile = 100.0 - lower_percentile

    baseline_mean = np.mean(boot_means, axis=0)
    baseline_std = np.std(boot_means, axis=0)
    baseline_ci_lower = np.percentile(boot_means, lower_percentile, axis=0)
    baseline_ci_upper = np.percentile(boot_means, upper_percentile, axis=0)
    return baseline_mean, baseline_std, baseline_ci_lower, baseline_ci_upper


def _compute_balanced_baseline(stationary_segments, running_segments, n_bootstraps=100, seed=42, ci_percentile=95):
    stationary_means = np.stack([seg.mean(axis=0) for seg in stationary_segments], axis=0)
    running_means = np.stack([seg.mean(axis=0) for seg in running_segments], axis=0)

    min_len = min(len(stationary_means), len(running_means))
    if min_len == 0:
        raise ValueError("Balanced passive baseline requires stationary and running windows.")

    rng = np.random.default_rng(seed)
    boot_means = []
    for _ in range(max(int(n_bootstraps), 1)):
        stat_idxs = rng.choice(len(stationary_means), size=min_len, replace=False)
        run_idxs = rng.choice(len(running_means), size=min_len, replace=False)
        sampled_means = np.concatenate([stationary_means[stat_idxs], running_means[run_idxs]], axis=0)
        boot_means.append(np.mean(sampled_means, axis=0))

    return _summarize_bootstrapped_means(boot_means, ci_percentile)


def _compute_pooled_baseline(segments, n_bootstraps=100, seed=42, ci_percentile=95):
    segment_means = np.stack([seg.mean(axis=0) for seg in segments], axis=0)
    if len(segment_means) == 0:
        raise ValueError("Pooled passive baseline requires at least one valid window.")

    rng = np.random.default_rng(seed)
    boot_means = []
    for _ in range(max(int(n_bootstraps), 1)):
        sample_idxs = rng.choice(len(segment_means), size=len(segment_means), replace=True)
        boot_means.append(np.mean(segment_means[sample_idxs], axis=0))

    return _summarize_bootstrapped_means(boot_means, ci_percentile)


sound_events, tl, cfg, paradigm, stimulus_catalog = _load_meta(snakemake.input[0])
if paradigm != "passive":
    raise ValueError("baseline_passive.py only supports passive sessions.")

with h5py.File(snakemake.input[1], "r") as f:
    lfp = np.array(f["lfp"]).T  # time x channels

with h5py.File(snakemake.input[2], "r") as f:
    artifact_mask = np.array(f["artifact_mask"]).astype(bool)

baseline_cfg = snakemake.config["lfp"]["baseline"]
passive_cfg = baseline_cfg.get("passive", {})

stationary_thresh = float(baseline_cfg["stationary_thresh"])
fs_lfp = float(snakemake.config["lfp"]["target_rate"])
post_offset_buffer_ms = float(passive_cfg.get("post_offset_buffer_ms", 100.0))
window_ms = float(passive_cfg.get("window_ms", 100.0))
n_bootstraps = int(passive_cfg.get("n_bootstraps", 100))
bootstrap_seed = int(passive_cfg.get("seed", 42))
ci_percentile = float(passive_cfg.get("ci_percentile", 95))

window_samples = max(int(round(window_ms * fs_lfp / 1000.0)), 1)
buffer_s = post_offset_buffer_ms / 1000.0
window_s = window_ms / 1000.0

n_samples = lfp.shape[0]
session_end_s = n_samples / fs_lfp

speed_times = tl[:, 0].astype(float)
speed_values = tl[:, 3].astype(float)
lfp_times = np.arange(n_samples, dtype=float) / fs_lfp
speed_interp = interp1d(
    speed_times,
    speed_values,
    kind="linear",
    bounds_error=False,
    fill_value=(float(speed_values[0]), float(speed_values[-1])),
)
speed_upsampled = speed_interp(lfp_times)

sound_events = sound_events[np.argsort(sound_events[:, 0])]
sound_ids = sound_events[:, 1].astype(np.int32)
has_explicit_silence = bool(np.any(sound_ids == 0))

positive_sound_ids = sorted({int(sound_id) for sound_id in sound_ids if int(sound_id) > 0})
missing_durations = [
    sound_id for sound_id in positive_sound_ids
    if _duration_sec(sound_id, stimulus_catalog) is None
]
if missing_durations:
    raise ValueError(
        "Passive baseline needs stimulus durations for all positive sound IDs. "
        f"Missing durations for: {missing_durations}"
    )

protected_mask = np.zeros(n_samples, dtype=bool)

for idx, row in enumerate(sound_events):
    onset_s = float(row[0])
    sound_id = int(row[1])

    if sound_id == 0:
        continue

    duration_s = _duration_sec(sound_id, stimulus_catalog)
    if duration_s is None:
        duration_s = max(_resolve_event_end_time(sound_events, idx, stimulus_catalog, session_end_s) - onset_s, 0.0)

    protected_end_s = onset_s + duration_s + buffer_s
    _mark_mask(protected_mask, onset_s, protected_end_s, fs_lfp)

stationary_segments = []
running_segments = []
valid_segments = []

accepted_windows = []
if has_explicit_silence:
    baseline_method = "passive_explicit_silence"
    rejection_counts = {
        "out_of_bounds": 0,
        "explicit_silence_too_short": 0,
        "artifact": 0,
    }

    last_stim_protected_end_s = -np.inf
    for idx, row in enumerate(sound_events):
        onset_s = float(row[0])
        sound_id = int(row[1])

        if sound_id > 0:
            duration_s = _duration_sec(sound_id, stimulus_catalog)
            if duration_s is None:
                duration_s = max(
                    _resolve_event_end_time(sound_events, idx, stimulus_catalog, session_end_s) - onset_s,
                    0.0,
                )
            last_stim_protected_end_s = max(last_stim_protected_end_s, onset_s + duration_s + buffer_s)
            continue

        silence_end_s = _resolve_event_end_time(sound_events, idx, stimulus_catalog, session_end_s)
        start_s = max(onset_s, last_stim_protected_end_s)
        end_s = start_s + window_s

        if end_s > session_end_s:
            rejection_counts["out_of_bounds"] += 1
            continue

        if end_s > silence_end_s:
            rejection_counts["explicit_silence_too_short"] += 1
            continue

        start_idx = _time_to_sample(start_s, fs_lfp, n_samples)
        end_idx = _time_to_sample(end_s, fs_lfp, n_samples)
        if end_idx <= start_idx or end_idx > n_samples:
            rejection_counts["out_of_bounds"] += 1
            continue

        if artifact_mask[start_idx:end_idx].any():
            rejection_counts["artifact"] += 1
            continue

        segment = lfp[start_idx:end_idx]
        speed_segment = speed_upsampled[start_idx:end_idx]
        if segment.shape[0] == 0:
            rejection_counts["out_of_bounds"] += 1
            continue

        speed_mean = float(np.nanmean(speed_segment))

        valid_segments.append(segment)
        if speed_mean < stationary_thresh:
            stationary_segments.append(segment)
            state_code = 0
        else:
            running_segments.append(segment)
            state_code = 1

        accepted_windows.append((sound_id, onset_s, start_s, end_s, speed_mean, state_code))
else:
    baseline_method = "passive_post_offset"
    rejection_counts = {
        "out_of_bounds": 0,
        "protected_overlap": 0,
        "artifact": 0,
    }

    for row in sound_events:
        onset_s = float(row[0])
        sound_id = int(row[1])
        if sound_id <= 0:
            continue

        duration_s = _duration_sec(sound_id, stimulus_catalog)
        start_s = onset_s + duration_s + buffer_s
        end_s = start_s + window_s
        if end_s > session_end_s:
            rejection_counts["out_of_bounds"] += 1
            continue

        start_idx = _time_to_sample(start_s, fs_lfp, n_samples)
        end_idx = _time_to_sample(end_s, fs_lfp, n_samples)
        if end_idx <= start_idx:
            rejection_counts["out_of_bounds"] += 1
            continue

        if end_idx > n_samples:
            rejection_counts["out_of_bounds"] += 1
            continue

        if protected_mask[start_idx:end_idx].any():
            rejection_counts["protected_overlap"] += 1
            continue

        if artifact_mask[start_idx:end_idx].any():
            rejection_counts["artifact"] += 1
            continue

        segment = lfp[start_idx:end_idx]
        speed_segment = speed_upsampled[start_idx:end_idx]
        if segment.shape[0] == 0:
            rejection_counts["out_of_bounds"] += 1
            continue

        speed_mean = float(np.nanmean(speed_segment))

        valid_segments.append(segment)
        if speed_mean < stationary_thresh:
            stationary_segments.append(segment)
            state_code = 0
        else:
            running_segments.append(segment)
            state_code = 1

        accepted_windows.append((sound_id, onset_s, start_s, end_s, speed_mean, state_code))

if stationary_segments and running_segments:
    baseline_mean, baseline_std, baseline_ci_lower, baseline_ci_upper = _compute_balanced_baseline(
        stationary_segments,
        running_segments,
        n_bootstraps=n_bootstraps,
        seed=bootstrap_seed,
        ci_percentile=ci_percentile,
    )
    state_mode = "balanced"
else:
    if not valid_segments:
        raise ValueError(
            "Passive baseline found no valid post-offset windows. "
            f"Rejections: {json.dumps(rejection_counts, sort_keys=True)}"
        )

    baseline_mean, baseline_std, baseline_ci_lower, baseline_ci_upper = _compute_pooled_baseline(
        valid_segments,
        n_bootstraps=n_bootstraps,
        seed=bootstrap_seed,
        ci_percentile=ci_percentile,
    )
    state_mode = "pooled"

with h5py.File(snakemake.output[0], "w") as f:
    base_mx = np.vstack([baseline_mean, baseline_std, baseline_ci_lower, baseline_ci_upper]).T
    ds = f.create_dataset("lfp_base", data=base_mx)
    ds.attrs["columns"] = "baseline_mean, baseline_std, baseline_ci_lower, baseline_ci_upper"

    f.attrs["baseline_method"] = baseline_method
    f.attrs["state_mode"] = state_mode
    f.attrs["post_offset_buffer_ms"] = post_offset_buffer_ms
    f.attrs["window_ms"] = window_ms
    f.attrs["uses_explicit_silence"] = int(has_explicit_silence)
    f.attrs["n_valid_segments"] = len(valid_segments)
    f.attrs["n_stationary_segments"] = len(stationary_segments)
    f.attrs["n_running_segments"] = len(running_segments)
    f.attrs["rejection_counts"] = json.dumps(rejection_counts, sort_keys=True)

    if accepted_windows:
        window_mx = np.array(accepted_windows, dtype=float)
        ds = f.create_dataset("baseline_windows", data=window_mx)
        ds.attrs["columns"] = "sound_id, stimulus_onset_s, baseline_start_s, baseline_end_s, mean_speed_mps, state_code"


fig, axes = plt.subplots(3, 1, figsize=(12, 8))

axes[0].errorbar(np.arange(len(baseline_mean)), baseline_mean, yerr=baseline_std, fmt=".")
axes[0].set_title("Passive baseline mean +/- std")

ci_width = baseline_ci_upper - baseline_ci_lower
axes[1].plot(ci_width, ".-")
axes[1].set_title("Passive baseline CI width per channel")

axes[2].axis("off")
summary_lines = [
    f"method: {baseline_method}",
    f"state_mode: {state_mode}",
    f"explicit silence events present: {has_explicit_silence}",
    f"post-offset buffer: {post_offset_buffer_ms:.1f} ms",
    f"window length: {window_ms:.1f} ms",
    f"valid windows: {len(valid_segments)}",
    f"stationary windows: {len(stationary_segments)}",
    f"running windows: {len(running_segments)}",
    f"rejections: {json.dumps(rejection_counts, sort_keys=True)}",
]
axes[2].text(0.01, 0.99, "\n".join(summary_lines), ha="left", va="top", family="monospace")

fig.tight_layout()
fig.savefig(snakemake.output[1])
plt.close(fig)
