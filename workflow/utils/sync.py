import numpy as np
from scipy import signal
from scipy.ndimage import uniform_filter1d
from scipy.signal import butter, detrend, sosfilt
from utils.neurosuite import XMLHero, DatHero
import os
import warnings

def _moving_average_power(x, n):
    n = max(1, int(n))
    p = x * x
    p_s = uniform_filter1d(p, size=n, mode="nearest")
    return np.sqrt(np.maximum(p_s, 0))


def _detect_onsets_hysteresis(env_s, hi, lo, min_on, min_off, min_gap):
    n = env_s.size
    onsets = []
    state_idle = True
    last_onset = -10**12

    i = 0
    while i < n:
        if state_idle:
            if i - last_onset >= min_gap and env_s[i] >= hi:
                end = min(n, i + min_on)
                if np.all(env_s[i:end] >= hi):
                    onsets.append(i)
                    last_onset = i
                    state_idle = False
                    i = end
                    continue
        else:
            if env_s[i] <= lo:
                end = min(n, i + min_off)
                if np.all(env_s[i:end] <= lo):
                    state_idle = True
                    i = end
                    continue
        i += 1

    return np.array(onsets, dtype=np.int64)


def find_all_sine_onsets_adc(
    t,
    x,
    fs,
    f_lo,
    f_hi,
    smooth_ms=1.0,
    k=6.0,
    hysteresis=0.6,
    min_on_ms=50.0,
    min_off_ms=3.0,
    min_gap_ms=50.0,
    baseline_frac=0.2,
    order=4,
):
    xd = detrend(x)

    sos = butter(order, [f_lo / (fs / 2), f_hi / (fs / 2)], btype="bandpass", output="sos")
    xf = sosfilt(sos, xd)

    n_smooth = max(1, int((smooth_ms / 1000) * fs))
    env_s = _moving_average_power(xf, n_smooth)

    n0 = max(10, int(baseline_frac * len(env_s)))
    base = env_s[:n0]
    med = np.median(base)
    mad = np.median(np.abs(base - med))
    sigma = 1.4826 * mad if mad > 0 else np.std(base)

    hi = med + k * sigma
    lo = med + (k * hysteresis) * sigma

    min_on = max(1, int((min_on_ms / 1000) * fs))
    min_off = max(1, int((min_off_ms / 1000) * fs))
    min_gap = max(0, int((min_gap_ms / 1000) * fs))

    onsets = _detect_onsets_hysteresis(env_s, hi, lo, min_on, min_off, min_gap)
    return onsets, t[onsets], env_s, (lo, hi)


def _format_timestamp_window(timestamps, center_idx, radius=3):
    start = max(0, center_idx - radius)
    stop = min(len(timestamps), center_idx + radius + 2)
    pairs = [f"{i}:{timestamps[i]:.12g}" for i in range(start, stop)]
    return "[" + ", ".join(pairs) + "]"


def _sample_numbers_diagnostic(timestamp_path):
    sample_numbers_path = os.path.join(os.path.dirname(timestamp_path), "sample_numbers.npy")
    if not os.path.exists(sample_numbers_path):
        return "sibling_sample_numbers: not found"

    sample_numbers = np.load(sample_numbers_path, mmap_mode="r")
    if len(sample_numbers) < 2:
        return (
            f"sibling_sample_numbers: {sample_numbers_path}\n"
            f"  count={len(sample_numbers)}; not enough values to assess continuity"
        )

    diffs = np.diff(sample_numbers)
    bad_idxs = np.where(diffs != 1)[0]
    examples = []
    for idx in bad_idxs[:8]:
        examples.append(
            f"{int(idx)}->{int(idx + 1)}: "
            f"{int(sample_numbers[idx])} -> {int(sample_numbers[idx + 1])} "
            f"(d={int(diffs[idx])})"
        )
    examples_text = ""
    if examples:
        examples_text = "\n  examples: " + "; ".join(examples)

    return (
        f"sibling_sample_numbers: {sample_numbers_path}\n"
        f"  count={len(sample_numbers)}, "
        f"min/max={int(np.min(sample_numbers))}/{int(np.max(sample_numbers))}, "
        f"non_contiguous_steps={len(bad_idxs)}"
        f"{examples_text}"
    )


def _load_contiguous_sample_numbers(timestamp_path, expected_len):
    sample_numbers_path = os.path.join(os.path.dirname(timestamp_path), "sample_numbers.npy")
    if not os.path.exists(sample_numbers_path):
        raise FileNotFoundError(
            f"Cannot rebuild timestamps because sample_numbers.npy is missing next to: "
            f"{timestamp_path}"
        )

    sample_numbers = np.load(sample_numbers_path, mmap_mode="r")
    if len(sample_numbers) != expected_len:
        raise ValueError(
            f"Cannot rebuild timestamps because sample_numbers length ({len(sample_numbers)}) "
            f"does not match timestamps length ({expected_len}): {sample_numbers_path}"
        )

    diffs = np.diff(sample_numbers)
    bad_idxs = np.where(diffs != 1)[0]
    if len(bad_idxs):
        examples = []
        for idx in bad_idxs[:8]:
            examples.append(
                f"{int(idx)}->{int(idx + 1)}: "
                f"{int(sample_numbers[idx])} -> {int(sample_numbers[idx + 1])} "
                f"(d={int(diffs[idx])})"
            )
        raise ValueError(
            f"Cannot rebuild timestamps because sample_numbers.npy is not contiguous: "
            f"{sample_numbers_path}\n"
            f"non_contiguous_steps={len(bad_idxs)}; examples: {'; '.join(examples)}"
        )

    return sample_numbers


def _repair_minor_timestamp_inversions(timestamps, path, label):
    timestamps = np.asarray(timestamps, dtype=float).copy()
    if len(timestamps) < 2:
        _validate_strictly_increasing_timestamps(timestamps, path, label)
        return timestamps

    diffs = np.diff(timestamps)
    bad_idxs = np.where(diffs <= 0)[0]
    if not len(bad_idxs):
        return timestamps

    positive_diffs = diffs[diffs > 0]
    median_dt = np.median(positive_diffs) if len(positive_diffs) else np.nan
    if not np.isfinite(median_dt) or median_dt <= 0:
        _validate_strictly_increasing_timestamps(timestamps, path, label)

    negative_value_count = int(np.sum(timestamps < 0))
    largest_backstep = float(np.max(-diffs[bad_idxs]))
    bad_fraction = len(bad_idxs) / max(1, len(diffs))

    # Open Ephys occasionally emits tiny one-sample timestamp inversions. Preserve
    # the OE clock in that case; only reject broad corruption such as -1 blocks.
    if (
        negative_value_count
        or bad_fraction > 1e-4
        or largest_backstep > 10 * float(median_dt)
    ):
        _validate_strictly_increasing_timestamps(timestamps, path, label)

    repaired = timestamps.copy()
    for idx in bad_idxs:
        start = int(idx + 1)
        repaired[start] = repaired[start - 1] + median_dt
        scan_idx = start + 1
        while scan_idx < len(repaired) and repaired[scan_idx] <= repaired[scan_idx - 1]:
            repaired[scan_idx] = repaired[scan_idx - 1] + median_dt
            scan_idx += 1

    warnings.warn(
        f"{label} timestamps.npy has {len(bad_idxs)} tiny non-increasing steps; "
        "preserving the Open Ephys clock and repairing those local inversions for "
        f"sync monotonicity. path: {path}",
        RuntimeWarning,
    )
    return repaired


def _load_adc_timestamps_for_sync(adc_ts_file, ephys_timestamps, s_rate):
    adc_timestamps = np.load(adc_ts_file)
    try:
        _validate_strictly_increasing_timestamps(
            adc_timestamps,
            adc_ts_file,
            "ADC",
        )
        return adc_timestamps, "timestamps.npy"
    except ValueError as exc:
        repaired = _repair_minor_timestamp_inversions(
            adc_timestamps,
            adc_ts_file,
            "ADC",
        )
        return repaired, "timestamps.npy repaired"


def _validate_strictly_increasing_timestamps(timestamps, path, label):
    timestamps = np.asarray(timestamps)

    problems = []
    finite_mask = np.isfinite(timestamps)
    nonfinite_idxs = np.where(~finite_mask)[0]
    if len(nonfinite_idxs):
        examples = ", ".join(str(int(i)) for i in nonfinite_idxs[:10])
        problems.append(
            f"non-finite values: count={len(nonfinite_idxs)}, first_indices=[{examples}]"
        )

    if len(timestamps) < 2:
        problems.append(f"not enough timestamps: count={len(timestamps)}")
        diffs = np.array([])
    else:
        diffs = np.diff(timestamps)
        bad_step_idxs = np.where(diffs <= 0)[0]
        if len(bad_step_idxs):
            neg_idxs = np.where(diffs < 0)[0]
            zero_idxs = np.where(diffs == 0)[0]
            examples = []
            for idx in bad_step_idxs[:8]:
                examples.append(
                    "idx "
                    f"{int(idx)}->{int(idx + 1)}: "
                    f"{timestamps[idx]:.12g} -> {timestamps[idx + 1]:.12g} "
                    f"(dt={diffs[idx]:.12g}); "
                    f"window={_format_timestamp_window(timestamps, int(idx))}"
                )
            neg_examples = []
            for idx in neg_idxs[:8]:
                neg_examples.append(
                    "idx "
                    f"{int(idx)}->{int(idx + 1)}: "
                    f"{timestamps[idx]:.12g} -> {timestamps[idx + 1]:.12g} "
                    f"(dt={diffs[idx]:.12g}); "
                    f"window={_format_timestamp_window(timestamps, int(idx))}"
                )
            neg_detail = ""
            if neg_examples:
                neg_detail = (
                    "\nfirst negative steps:\n  - " + "\n  - ".join(neg_examples)
                )
            problems.append(
                "non-increasing steps: "
                f"count={len(bad_step_idxs)}, "
                f"negative={len(neg_idxs)}, zero={len(zero_idxs)}, "
                "examples:\n  - " + "\n  - ".join(examples)
                + neg_detail
            )

    if problems:
        neg_value_count = int(np.sum(timestamps < 0)) if len(timestamps) else 0
        min_value = np.nanmin(timestamps) if len(timestamps) else np.nan
        max_value = np.nanmax(timestamps) if len(timestamps) else np.nan
        median_dt = np.nanmedian(diffs) if len(diffs) else np.nan
        raise ValueError(
            f"Invalid {label} timestamp vector for sound/ephys sync.\n"
            f"path: {path}\n"
            f"count: {len(timestamps)}\n"
            f"min/max: {min_value:.12g} / {max_value:.12g}\n"
            f"negative_value_count: {neg_value_count}\n"
            f"median_dt: {median_dt:.12g}\n"
            f"{_sample_numbers_diagnostic(path)}\n"
            + "\n".join(problems)
        )


def _merge_onsets_into_events(
    events,
    onsets_s,
    max_missed_gap_s=0.1,
    anomaly_gap_s=-10.0,
):
    events = events.copy()
    detectable_mask = events[:, 1] > 0
    detectable_indices = np.where(detectable_mask)[0]
    detectable_events = events[detectable_mask].copy()

    merged_events = events.copy()
    excluded_detectable_indices = []
    extra_excluded_indices = []
    count_missed = 0

    onsets_s = np.sort(np.asarray(onsets_s))

    for t_i_ana, t_ana in enumerate(onsets_s):
        expected_idx = t_i_ana - count_missed
        if expected_idx >= len(detectable_events):
            break

        gap = t_ana - detectable_events[expected_idx, 0]
        original_idx = detectable_indices[expected_idx]

        if gap > max_missed_gap_s:
            excluded_detectable_indices.append(original_idx)
            count_missed -= 1
            continue

        if gap < anomaly_gap_s:
            excluded_detectable_indices.append(original_idx)
            count_missed += 1
            continue

        if gap < 0:
            detectable_events[expected_idx:, 0] += gap
            merged_events[original_idx:, 0] += gap

        merged_events[original_idx, 0] = t_ana

        if original_idx > 0 and merged_events[original_idx, 0] < merged_events[original_idx - 1, 0]:
            prev_idx = original_idx - 1
            while prev_idx >= 0 and merged_events[original_idx, 0] < merged_events[prev_idx, 0]:
                if merged_events[prev_idx, 1] <= 0:
                    extra_excluded_indices.append(prev_idx)
                    prev_idx -= 1
                    continue
                excluded_detectable_indices.append(original_idx)
                break

    excluded_indices = sorted(set(excluded_detectable_indices + extra_excluded_indices))
    if excluded_indices:
        keep_mask = np.ones(len(merged_events), dtype=bool)
        keep_mask[excluded_indices] = False
        merged_events = merged_events[keep_mask]

    return merged_events



def infer_adc_ch_no(adc_file, adc_ts_file, dtype=np.int16):
    """
    Infer number of ADC channels from file size and timestamps length.
    Assumes timestamps.npy has one entry per sample frame.
    """
    adc_timestamps = np.load(adc_ts_file)
    if len(adc_timestamps) == 0:
        raise ValueError(f"ADC timestamps file is empty: {adc_ts_file}")

    n_values = os.path.getsize(adc_file) // np.dtype(dtype).itemsize
    if n_values % len(adc_timestamps) != 0:
        raise ValueError(
            f"Cannot infer ADC channel count: total values={n_values} not divisible "
            f"by number of timestamps={len(adc_timestamps)} for {adc_file}"
        )

    ch_no = n_values // len(adc_timestamps)
    if ch_no <= 0:
        raise ValueError(f"Inferred invalid ADC channel count {ch_no} for {adc_file}")
    return int(ch_no)


def _load_adc_channel(adc_file, adc_ts_file, channel, ch_no=None, dtype=np.int16):
    if ch_no is None:
        ch_no = infer_adc_ch_no(adc_file, adc_ts_file, dtype=dtype)

    if channel >= ch_no:
        raise ValueError(
            f"Requested ADC channel {channel}, but channel count is only {ch_no} "
            f"for file {adc_file}"
        )

    raw = np.memmap(adc_file, dtype=dtype, mode="r")
    if raw.size % ch_no != 0:
        raise ValueError(
            f"ADC file size ({raw.size} values) is not divisible by channel count {ch_no} "
            f"for file {adc_file}"
        )

    n_frames = raw.size // ch_no
    adc_timestamps = np.load(adc_ts_file, mmap_mode="r")
    if n_frames != len(adc_timestamps):
        raise ValueError(
            f"ADC frame count mismatch: dat implies {n_frames} frames, "
            f"timestamps has {len(adc_timestamps)} entries"
        )

    data = raw.reshape(n_frames, ch_no)
    return np.asarray(data[:, channel], dtype=float), int(ch_no), int(n_frames)


def get_sound_events_from_ADC(
    adc_file,
    adc_ts_file,
    sounds_file,
    events_file,
    ephys_ts_file,
    channel=11,
    event_th=300,
    s_rate=None,
    ch_no=None,
    dtype=np.int16,
):
    """
    Detect sound pulse periods from ADC channel and align logged sounds to ephys time.

    Parameters
    ----------
    adc_file : str
        Path to ADC .dat file
    adc_ts_file : str
        Path to ADC timestamps.npy
    sounds_file : str
        Logged sounds CSV
    events_file : str
        Logged events CSV
    ephys_ts_file : str
        Timestamps of a reference probe stream
    channel : int
        ADC channel carrying sound pulses
    event_th : float
        Threshold after smoothing
    s_rate : float or None
        ADC sample rate. If None, inferred from adc timestamps.
    ch_no : int or None
        Number of ADC channels. If None, inferred from file size and timestamps.
    dtype : numpy dtype
        ADC binary dtype
    """

    # --- load metadata first ---
    ephys_timestamps = np.load(ephys_ts_file)
    events_exp = np.loadtxt(events_file, skiprows=1, delimiter=",")
    events_csv = np.loadtxt(sounds_file, skiprows=1, delimiter=",")

    if events_csv.ndim == 1:
        events_csv = events_csv[None, :]
    if events_exp.ndim == 1:
        events_exp = events_exp[None, :]

    # normalize logged sounds to session start
    events_csv[:, 0] = events_csv[:, 0] - events_exp[0][0]

    # --- infer channel count if needed ---
    if ch_no is None:
        ch_no = infer_adc_ch_no(adc_file, adc_ts_file, dtype=dtype)

    if channel >= ch_no:
        raise ValueError(
            f"Requested ADC channel {channel}, but inferred channel count is only {ch_no} "
            f"for file {adc_file}"
        )

    # --- infer sample rate if needed ---
    if s_rate is None:
        adc_timestamps = np.load(adc_ts_file)
        _validate_strictly_increasing_timestamps(
            adc_timestamps,
            adc_ts_file,
            "ADC",
        )
        if len(adc_timestamps) < 2:
            raise ValueError("Cannot infer ADC sample rate from fewer than 2 timestamps")
        dt = np.median(np.diff(adc_timestamps))
        if dt <= 0:
            raise ValueError(f"Invalid ADC timestamp spacing: median dt={dt}")
        s_rate = 1.0 / dt

    adc_timestamps, adc_clock_source = _load_adc_timestamps_for_sync(
        adc_ts_file,
        ephys_timestamps,
        s_rate,
    )

    # --- load ADC data robustly ---
    channel_data, ch_no, n_frames = _load_adc_channel(
        adc_file,
        adc_ts_file,
        channel,
        ch_no=ch_no,
        dtype=dtype,
    )

    # --- smoothing and thresholding ---
    kernel_width = max(3, int(round(s_rate / 100)))  # ~10 ms window
    kernel = signal.windows.gaussian(kernel_width, std=kernel_width / 7.2)
    kernel /= kernel.sum()

    data_smooth = np.convolve(
        np.abs(channel_data - channel_data.mean()),
        kernel,
        mode="same",
    )

    idxs_high = np.where(data_smooth > event_th)[0]

    if len(idxs_high) == 0:
        # return empty detected periods + unchanged logged events
        return np.zeros((0, 2)), events_csv.copy()

    # --- detect contiguous high periods ---
    periods = []
    gap_breaks = np.where(np.diff(idxs_high) > 3)[0]

    start_idx = idxs_high[0]
    for break_idx in gap_breaks:
        end_idx = idxs_high[break_idx]
        periods.append((start_idx, end_idx))
        start_idx = idxs_high[break_idx + 1]

    # append final period
    periods.append((start_idx, idxs_high[-1]))
    periods = np.asarray(periods, dtype=int)

    # --- convert sample indices to times relative to ephys start ---
    p_begs = adc_timestamps[periods[:, 0]] - ephys_timestamps[0]
    p_ends = adc_timestamps[periods[:, 1]] - ephys_timestamps[0]
    period_times = np.column_stack([p_begs, p_ends])

    # --- align logged sound events to detected ADC pulse times ---
    adc_to_ephys_shift = adc_timestamps[0] - ephys_timestamps[0]
    shift = events_csv[0][0] - period_times[0][0] + adc_to_ephys_shift
    events_synced = events_csv.copy()
    events_synced[:, 0] = events_synced[:, 0] - shift

    # linear drift correction
    t_last_ADC = period_times[-1][0]
    nearest_idx = np.abs(t_last_ADC - events_synced[:, 0]).argmin()
    t_last_logged = events_synced[nearest_idx][0]
    drift = t_last_ADC - t_last_logged
    events_synced[:, 0] = events_synced[:, 0] + np.arange(len(events_synced)) * drift / len(events_synced)

    return period_times, events_synced

def refine_sound_events_from_ADC(
    period_times,
    events_synced,
    adc_file,
    adc_ts_file,
    ephys_ts_file,
    channel=11,
    s_rate=30300,
    ch_no=None,
    dtype=np.int16,
    f_lo=600.0,
    f_hi=1400.0,
    smooth_ms=1.0,
    k=6.0,
    hysteresis=0.6,
    min_on_ms=50.0,
    min_off_ms=3.0,
    min_gap_ms=50.0,
    baseline_frac=0.2,
    order=4,
    max_missed_gap_s=0.1,
    anomaly_gap_s=-10.0,
):
    channel_data, ch_no, n_frames = _load_adc_channel(
        adc_file,
        adc_ts_file,
        channel,
        ch_no=ch_no,
        dtype=dtype,
    )

    ephys_timestamps = np.load(ephys_ts_file)
    adc_timestamps, adc_clock_source = _load_adc_timestamps_for_sync(
        adc_ts_file,
        ephys_timestamps,
        s_rate,
    )
    t = adc_timestamps - ephys_timestamps[0]

    _, onsets_s, _, _ = find_all_sine_onsets_adc(
        t=t,
        x=channel_data,
        fs=s_rate,
        f_lo=f_lo,
        f_hi=f_hi,
        smooth_ms=smooth_ms,
        k=k,
        hysteresis=hysteresis,
        min_on_ms=min_on_ms,
        min_off_ms=min_off_ms,
        min_gap_ms=min_gap_ms,
        baseline_frac=baseline_frac,
        order=order,
    )

    events_refined = _merge_onsets_into_events(
        events_synced,
        onsets_s,
        max_missed_gap_s=max_missed_gap_s,
        anomaly_gap_s=anomaly_gap_s,
    )

    return period_times, events_refined

# def get_sound_events_from_ADC(adc_file, adc_ts_file, sounds_file, events_file, ephys_ts_file, channel=11, event_th=300, s_rate=30300):
#     # read ADC channel with sound pulses
#     dh = DatHero(adc_file, s_rate=s_rate, ch_no=12)
#     channel_data = dh.get_single_channel(channel)

#     # smoothing and thresholding
#     kernel_width = int(s_rate / 100)  # need to test if good enough for high frequencies
#     kernel = signal.windows.gaussian(kernel_width, std=(kernel_width) / 7.2)

#     data_smooth = np.convolve(np.abs(channel_data - channel_data.mean()), kernel, 'same') / kernel.sum()

#     # TODO: make threshold dependent on noise levels between events
#     idxs_high = np.where(data_smooth > event_th)[0]  # indices where sound was ON
#     if len(idxs_high) == 0:
#         return np.zeros((0,2)), np.zeros((0, events_csv.shape[1]))
    
#     # detect sound events
#     periods = []
#     idxs_diff = np.diff(idxs_high)
#     period_idxs = np.where((idxs_diff > 3))[0]

#     for i, idx in enumerate(period_idxs):
#         if i == 0:
#             pair = (idxs_high[0], idxs_high[idx])
#         else:
#             pair = (idxs_high[period_idxs[i - 1] + 1], idxs_high[idx])
#         periods.append(pair)
#     periods = np.array(periods)  # these are all pairs of sample indices

#     # convert samples indices into actual times in seconds, with the zero as the start of ephys recording
#     # start of ephys recording is NOT the same as the start of ADC recording
#     adc_timestamps = np.load(adc_ts_file)
#     ephys_timestamps = np.load(ephys_ts_file)
#     p_begs = adc_timestamps[periods[:, 0]] - ephys_timestamps[0]  # start of the period in seconds
#     p_ends = adc_timestamps[periods[:, 1]] - ephys_timestamps[0]  # end of the period in seconds
#     period_times = np.column_stack([p_begs, p_ends])

#     events_exp = np.loadtxt(events_file, skiprows=1, delimiter=',')
#     events_csv = np.loadtxt(sounds_file, skiprows=1, delimiter=',')
#     events_csv[:, 0] = events_csv[:, 0] - events_exp[0][0]

#     # first shift by the delay of the ephys start
#     adc_to_ephys_shift = adc_timestamps[0] - ephys_timestamps[0]
#     shift = events_csv[0][0] - period_times[0][0] + adc_to_ephys_shift
#     events_csv[:, 0] = events_csv[:, 0] - shift

#     # next linearly correct the drift by finding a time diff between
#     # the last ADC pulse and the last logged one, assuming they are still the closest events
#     t_last_ADC = period_times[-1][0]
#     t_last_sev = events_csv[np.abs(t_last_ADC - events_csv[:, 0]).argmin()][0]
#     drift = t_last_ADC - t_last_sev
#     events_csv[:, 0] = events_csv[:, 0] + np.arange(len(events_csv)) * drift/len(events_csv)

#     return period_times, events_csv


def get_sound_events_from_openephys(dat_file, xml_file, sounds_file, events_file, channel, event_th=200, ipi=0.25, bias_correction=0.004):
    """
    returns: 
     - events_detected - sound events (t_start, t_end) detected from ephys
     - events_synced   - sound events (t_start, type) synced with the logger
    """
    # read ephys events channel
    xml_hero = XMLHero(xml_file)
    s_rate   = xml_hero.get_sampling_rate()
    ch_count = xml_hero.get_channel_count()

    dat_hero = DatHero(dat_file, s_rate=s_rate, ch_no=ch_count)

    data = dat_hero.get_single_channel(channel_no=channel)

    # smoothing and thresholding
    kernel_width = s_rate / 100  # need to test if good enough for high frequencies
    kernel = signal.gaussian(kernel_width, std=(kernel_width) / 7.2)

    data_smooth = np.convolve(np.abs(data), kernel, 'same') / kernel.sum()

    # TODO: make threshold dependent on noise levels between events
    idxs_high = np.where(data_smooth > event_th)[0]
    
    # detect sound events
    periods = []

    idxs_diff = np.diff(idxs_high)
    period_idxs = np.where((idxs_diff > 3))[0]

    for i, idx in enumerate(period_idxs):
        if i == 0:
            pair = (idxs_high[0], idxs_high[idx])
        else:
            pair = (idxs_high[period_idxs[i - 1] + 1], idxs_high[idx])
        periods.append(pair)

    # add the last period (ignoring events that were still happening at the end of recording)
    if not idxs_high[-1] == len(data_smooth):
        periods.append((idxs_high[period_idxs[-1] + 1], idxs_high[-1]))

    periods = np.array(periods)/float(s_rate)

    durations = np.diff(periods, axis=1)
    durations = durations.T[0]
    
    # sync events
    max_drift = 0.1  # maximum drift between logged and ephys events
    max_pulse = 0.22 # maximum pulse duration
    min_noise = 1    # at least one second for noise

    events_exp = np.loadtxt(events_file, skiprows=1, delimiter=',')
    events_csv = np.loadtxt(sounds_file, skiprows=1, delimiter=',')
    events_csv[:, 0] = events_csv[:, 0] - events_exp[0][0]

    events_synced = events_csv.copy()

    pulses = periods[durations < max_pulse]
    noises = periods[durations > min_noise]

    # update trial (BGR / TGT) pulse times
    for pulse in pulses:
        idx_ev = np.abs(pulse[0] - events_csv[:, 0]).argmin()
        if np.abs(pulse[0] - events_csv[idx_ev][0]) < max_drift:
            events_synced[idx_ev][0] = pulse[0] - bias_correction # correct for 4ms bias in the pulse detection

    # update noise periods
    for noise in noises:
        idx_ev = np.abs(noise[0] - events_csv[:, 0]).argmin()
        curr_drift = noise[0] - events_csv[idx_ev][0]
        if curr_drift < max_drift:
            idxs = np.where((events_csv[:, 0] > noise[0] - curr_drift - 0.01) & (events_csv[:, 0] < noise[1]))[0]
            for idx in idxs:
                events_synced[idx, 0] += curr_drift
                events_synced[idx, 1] = -1

    # update no stimulus periods
    curr_drift = None
    for i, event in enumerate(events_synced[1:]):
        if events_synced[i][1] == 0 and not events_synced[i-1][1] == 0:
            curr_drift = events_synced[i-1][0] - (events_synced[i][0] - ipi)
        if events_synced[i][1] == 0:
            events_synced[i][0] += curr_drift
            
    return periods, events_synced
