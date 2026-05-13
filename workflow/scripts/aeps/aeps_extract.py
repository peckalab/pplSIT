import os
import sys
import json
import h5py
import numpy as np

# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.lfp import notch_50hz, highpass_lfp, compute_aep_snr


def _is_legacy_aeps_format(aeps_cfg: dict) -> bool:
    """
    Legacy format:
    {
        "AEPs": {
            "A1": [3, 5, 7],
            "PPC": [20, 21]
        }
    }

    New format:
    {
        "AEPs": {
            "ProbeA": {
                "A1": {...},
                "PPC": {...}
            }
        }
    }
    """
    if not isinstance(aeps_cfg, dict) or len(aeps_cfg) == 0:
        return False

    first_value = next(iter(aeps_cfg.values()))
    return isinstance(first_value, (list, int))


def _resolve_current_stream_name(aeps_cfg: dict, snakemake) -> str | None:
    """
    Resolve the current stream name for new-format AEP configs.
    Returns None for legacy format.
    """
    if _is_legacy_aeps_format(aeps_cfg):
        return None

    # Preferred: Snakemake wildcard
    if hasattr(snakemake, "wildcards") and hasattr(snakemake.wildcards, "stream"):
        stream_name = str(snakemake.wildcards.stream)
        if stream_name in aeps_cfg:
            return stream_name

    # Fallback: only one stream present in the config
    if len(aeps_cfg) == 1:
        return next(iter(aeps_cfg.keys()))

    raise ValueError(
        "Could not resolve current stream name for new-format AEPs config. "
        "Expected snakemake.wildcards.stream or exactly one stream in manual.json."
    )


def _resolve_area_specs(aeps_cfg: dict, snakemake) -> tuple[dict, str]:
    """
    Returns:
        area_specs: dict mapping area -> spec
        stream_name: resolved stream name, or 'legacy'
    """
    if _is_legacy_aeps_format(aeps_cfg):
        return aeps_cfg, "legacy"

    stream_name = _resolve_current_stream_name(aeps_cfg, snakemake)
    if stream_name not in aeps_cfg:
        raise KeyError(f"Stream '{stream_name}' not found in manual.json['AEPs'].")

    area_specs = aeps_cfg[stream_name]
    if not isinstance(area_specs, dict) or len(area_specs) == 0:
        raise ValueError(f"No area specs found for stream '{stream_name}'.")

    return area_specs, stream_name


def _normalize_legacy_area_spec(spec):
    """
    Convert legacy per-area spec into a unified new-style spec.
    Examples:
        [3, 5, 7] -> {"channels": [3, 5, 7], "select": {"method": "all"}}
        12        -> {"channels": [12],       "select": {"method": "all"}}
    """
    if isinstance(spec, int):
        channels = [spec]
    elif isinstance(spec, list):
        channels = spec
    else:
        raise ValueError(
            f"Legacy AEP area spec must be int or list[int], got: {type(spec)}"
        )

    return {
        "channels": channels,
        "select": {"method": "all"},
    }


def _resolve_candidate_channels(channels_spec, n_channels_total: int) -> np.ndarray:
    """
    Supported:
      - "all"
      - [ch1, ch2, ch3, ...] explicit channel list
      - single int
    """
    if channels_spec == "all":
        candidate_channels = np.arange(n_channels_total, dtype=int)

    elif isinstance(channels_spec, int):
        candidate_channels = np.array([channels_spec], dtype=int)

    elif isinstance(channels_spec, list):
        candidate_channels = np.asarray(channels_spec, dtype=int)

    else:
        raise ValueError(
            "Unsupported 'channels' spec. Use 'all', int, or list[int]."
        )

    if candidate_channels.size == 0:
        raise ValueError("Candidate channel pool is empty.")

    if np.any(candidate_channels < 0) or np.any(candidate_channels >= n_channels_total):
        raise ValueError(
            f"Candidate channels out of bounds. Valid range: [0, {n_channels_total - 1}]"
        )

    return np.unique(candidate_channels)


def _select_channels(
    lfp_clean: np.ndarray,
    stim_times: np.ndarray,
    fs: float,
    candidate_channels: np.ndarray,
    select_spec: dict,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Returns
    -------
    selected_channels : np.ndarray
    channel_scores : np.ndarray
        Score for each candidate channel, same order as candidate_channels.
    channel_evoked_ptp : np.ndarray
    channel_sustained_std : np.ndarray
    """
    method = select_spec.get("method", "all")

    channel_scores = np.full(candidate_channels.shape, np.nan, dtype=float)
    channel_evoked_ptp = np.full(candidate_channels.shape, np.nan, dtype=float)
    channel_sustained_std = np.full(candidate_channels.shape, np.nan, dtype=float)

    if method == "all":
        selected_channels = candidate_channels.copy()

    elif method == "top_aep_snr":
        n_select = int(select_spec.get("n", 10))
        if n_select <= 0:
            raise ValueError("'select.n' must be a positive integer for top_aep_snr.")

        for i, ch in enumerate(candidate_channels):
            res = compute_aep_snr(
                signal_1d=lfp_clean[:, ch],
                stim_times=stim_times,
                fs=fs,
            )
            channel_scores[i] = res["score"]
            channel_evoked_ptp[i] = res["evoked_ptp"]
            channel_sustained_std[i] = res["sustained_std"]

        finite_mask = np.isfinite(channel_scores)
        if not np.any(finite_mask):
            raise ValueError(
                "No finite AEP-SNR scores were computed for candidate channel pool."
            )

        finite_candidates = candidate_channels[finite_mask]
        finite_scores = channel_scores[finite_mask]

        order = np.argsort(finite_scores)[::-1]  # descending
        n_take = min(n_select, finite_candidates.size)
        selected_channels = finite_candidates[order[:n_take]]

    else:
        raise ValueError(
            f"Unknown selection method '{method}'. Supported: 'all', 'top_aep_snr'."
        )

    if selected_channels.size == 0:
        raise ValueError("No channels selected for AEP extraction.")

    return selected_channels, channel_scores, channel_evoked_ptp, channel_sustained_std


def _extract_area_aeps(
    lfp_clean: np.ndarray,
    lfp_dirty: np.ndarray,
    stim_times: np.ndarray,
    selected_channels: np.ndarray,
    baseline_stds: np.ndarray,
    amplitude_weights: np.ndarray,
    post_samp: int,
) -> dict:
    """
    Extract AEPs for one area.

    Notes
    -----
    - Keeps your current design: post-stimulus only, no pre-stimulus extraction.
    - Weighted average uses per-area normalized weights.
    - Ensures number of extracted AEPs exactly matches number of sound_events:
      if a segment would overflow recording bounds, the missing tail is padded with zeros.
    """
    if selected_channels.size == 0:
        raise ValueError("selected_channels is empty.")

    # Renormalize weights within selected channels only
    selected_amp_weights = np.asarray(amplitude_weights[selected_channels], dtype=float)
    selected_amp_weights = np.nan_to_num(
        selected_amp_weights, nan=0.0, posinf=0.0, neginf=0.0
    )

    if np.all(selected_amp_weights <= 0):
        local_weights = np.ones(selected_channels.size, dtype=float) / selected_channels.size
    else:
        local_weights = selected_amp_weights / (np.sum(selected_amp_weights) + 1e-12)

    selected_baseline_stds = np.asarray(baseline_stds[selected_channels], dtype=float)
    selected_baseline_stds = np.where(
        np.isfinite(selected_baseline_stds) & (selected_baseline_stds > 0),
        selected_baseline_stds,
        np.nan,
    )

    avg_across_channels = []
    weighted_avg = []
    zscored_avg = []
    raw_AEPs = []

    n_time = lfp_clean.shape[0]
    n_sel = selected_channels.size

    for t in stim_times:
        t = int(round(t))
        start = t
        end = t + post_samp

        # Always allocate a full-length segment so output matches sound_events
        seg = np.zeros((post_samp, n_sel), dtype=lfp_clean.dtype)

        # Copy valid overlap from recording into the segment
        src_start = max(start, 0)
        src_end = min(end, n_time)

        if src_end > src_start:
            dst_start = src_start - start
            dst_end = dst_start + (src_end - src_start)
            seg[dst_start:dst_end, :] = lfp_clean[src_start:src_end, selected_channels]

        # Keep current behavior: no baseline subtraction
        seg_baselined = seg

        avg_across_channels.append(np.nanmean(seg_baselined, axis=1))
        weighted_avg.append(np.nansum(seg_baselined * local_weights[None, :], axis=1))

        seg_z = seg_baselined / selected_baseline_stds[None, :]
        zscored_avg.append(np.nanmean(seg_z, axis=1))

        raw_AEPs.append(seg_baselined)

    avg_across_channels = np.stack(avg_across_channels, axis=0)
    weighted_avg = np.stack(weighted_avg, axis=0)
    zscored_avg = np.stack(zscored_avg, axis=0)
    raw_AEPs = np.stack(raw_AEPs, axis=0)

    return {
        "avg_across_channels": avg_across_channels,
        "weighted_avg": weighted_avg,
        "zscored_avg": zscored_avg,
        "raw_AEPs": raw_AEPs,  # pulse, samples, channels
        "lfp_dirty": np.mean(lfp_dirty[:, selected_channels], axis=1),
        "lfp_clean": np.nanmean(lfp_clean[:, selected_channels], axis=1),
        "local_weights": local_weights,
    }


def _write_area_group(
    h5_group,
    area_name: str,
    area_result: dict,
    candidate_channels: np.ndarray,
    selected_channels: np.ndarray,
    channel_scores: np.ndarray,
    channel_evoked_ptp: np.ndarray,
    channel_sustained_std: np.ndarray,
    selection_method: str,
    stream_name: str,
):
    if area_name in h5_group:
        del h5_group[area_name]

    g = h5_group.create_group(area_name)

    g.create_dataset("avg_across_channels", data=area_result["avg_across_channels"])
    g.create_dataset("weighted_avg", data=area_result["weighted_avg"])
    g.create_dataset("zscored_avg", data=area_result["zscored_avg"])
    g.create_dataset("raw_AEPs", data=area_result["raw_AEPs"])

    g.create_dataset("lfp_dirty", data=area_result["lfp_dirty"])
    g.create_dataset("lfp_clean", data=area_result["lfp_clean"])

    g.create_dataset("candidate_channels", data=candidate_channels.astype(int))
    g.create_dataset("selected_channels", data=selected_channels.astype(int))
    g.create_dataset("selection_metric", data=np.asarray(channel_scores, dtype=float))
    g.create_dataset("selection_evoked_ptp", data=np.asarray(channel_evoked_ptp, dtype=float))
    g.create_dataset("selection_sustained_std", data=np.asarray(channel_sustained_std, dtype=float))
    g.create_dataset("selection_weights", data=np.asarray(area_result["local_weights"], dtype=float))

    g.attrs["selection_method"] = selection_method
    g.attrs["stream_name"] = stream_name


# -----------------------------------------------------------------------------
# Load inputs
# -----------------------------------------------------------------------------

with open(snakemake.input[0]) as json_file:
    manual_cfg = json.load(json_file)

if "AEPs" not in manual_cfg:
    raise KeyError("manual.json must contain the key 'AEPs'.")

aeps_cfg = manual_cfg["AEPs"]

with h5py.File(snakemake.input[1], "r") as f:
    sound_events = np.array(f["processed"]["sound_events"])

with h5py.File(snakemake.input[2], "r") as f:
    lfp = np.array(f["lfp"]).T  # time x channels

with h5py.File(snakemake.input[3], "r") as f:
    baselines = np.array(f["lfp_base"])
    baseline_per_channel = baselines[:, 0]  # kept for compatibility/future use
    baseline_stds = baselines[:, 1]

with h5py.File(snakemake.input[4], "r") as f:
    artifact_mask = np.array(f["artifact_mask"]).astype(bool)

with h5py.File(snakemake.input[5], "r") as f:
    metrics = np.array(f["aeps_lfp_metrics"])  # e.g. snrs, amplitudes, pre-stim baselines
    amplitude_weights = metrics[:, 1]

fs = snakemake.config["lfp"]["target_rate"]  # Hz
pre_ms = snakemake.config["lfp"]["aep_pre"] * 1000  # kept, but not used in extraction by design
post_ms = snakemake.config["lfp"]["aep_dur"] * 1000
signal_window_ms = snakemake.config["lfp"]["aep_sig"] * 1000  # kept for compatibility

stim_times = sound_events[:, 0] * fs

n_ch = lfp.shape[1]
pre_samp = int(pre_ms * fs / 1000)  # intentionally unused in extraction
post_samp = int(post_ms * fs / 1000)
signal_win = int(signal_window_ms * fs / 1000)  # intentionally unused
win_len = pre_samp + post_samp  # intentionally unused

if baseline_stds.shape[0] != n_ch:
    raise ValueError("baseline_stds length does not match number of LFP channels.")
if amplitude_weights.shape[0] != n_ch:
    raise ValueError("amplitude_weights length does not match number of LFP channels.")

# -----------------------------------------------------------------------------
# Clean LFP
# -----------------------------------------------------------------------------

lfp = notch_50hz(lfp, fs=fs, freq=50.0, q=30.0)
lfp_dirty = highpass_lfp(lfp, fs=fs, cutoff=2.0, order=4)

lfp_std = np.nanstd(lfp_dirty)
lfp_dirty = np.clip(lfp_dirty, -5 * lfp_std, 5 * lfp_std)

lfp_clean = lfp_dirty.copy()
lfp_clean[artifact_mask, :] = np.nan

# -----------------------------------------------------------------------------
# Resolve AEP area specs
# -----------------------------------------------------------------------------

area_specs, stream_name = _resolve_area_specs(aeps_cfg, snakemake)

# -----------------------------------------------------------------------------
# Write output
# -----------------------------------------------------------------------------

with h5py.File(snakemake.output[0], "w") as f_out:
    root_group = f_out

    for area, raw_spec in area_specs.items():
        # Convert legacy area specs to unified style
        if _is_legacy_aeps_format({area: raw_spec}):
            spec = _normalize_legacy_area_spec(raw_spec)
        else:
            spec = raw_spec

        if not isinstance(spec, dict):
            raise ValueError(f"Area spec for '{area}' must be a dict after normalization.")

        if "channels" not in spec:
            raise KeyError(f"Area '{area}' spec must contain 'channels'.")

        candidate_channels = _resolve_candidate_channels(
            spec["channels"], n_channels_total=n_ch
        )

        select_spec = spec.get("select", {"method": "all"})
        selection_method = select_spec.get("method", "all")

        selected_channels, channel_scores, channel_evoked_ptp, channel_sustained_std = _select_channels(
            lfp_clean=lfp_clean,
            stim_times=stim_times,
            fs=fs,
            candidate_channels=candidate_channels,
            select_spec=select_spec,
        )

        area_result = _extract_area_aeps(
            lfp_clean=lfp_clean,
            lfp_dirty=lfp_dirty,
            stim_times=stim_times,
            selected_channels=selected_channels,
            baseline_stds=baseline_stds,
            amplitude_weights=amplitude_weights,
            post_samp=post_samp,
        )

        _write_area_group(
            h5_group=root_group,
            area_name=area,
            area_result=area_result,
            candidate_channels=candidate_channels,
            selected_channels=selected_channels,
            channel_scores=channel_scores,
            channel_evoked_ptp=channel_evoked_ptp,
            channel_sustained_std=channel_sustained_std,
            selection_method=selection_method,
            stream_name=stream_name,
        )