import numpy as np
from scipy import signal
from scipy.ndimage import uniform_filter1d
from scipy.signal import butter, detrend, sosfilt
from utils.neurosuite import XMLHero, DatHero


def _moving_average_power(x, n):
    n = max(1, int(n))
    p = x * x
    p_s = uniform_filter1d(p, size=n, mode="nearest")
    return np.sqrt(p_s)


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


def refine_sound_events_from_ADC(
    period_times,
    events_synced,
    adc_file,
    adc_ts_file,
    ephys_ts_file,
    channel=11,
    s_rate=30300,
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
    dh = DatHero(adc_file, s_rate=s_rate, ch_no=12)
    channel_data = dh.get_single_channel(channel)

    adc_timestamps = np.load(adc_ts_file)
    ephys_timestamps = np.load(ephys_ts_file)
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


def get_sound_events_from_ADC(adc_file, adc_ts_file, sounds_file, events_file, ephys_ts_file, channel=11, event_th=300, s_rate=30300):
    # read ADC channel with sound pulses
    dh = DatHero(adc_file, s_rate=s_rate, ch_no=12)
    channel_data = dh.get_single_channel(channel)

    # smoothing and thresholding
    kernel_width = s_rate / 100  # need to test if good enough for high frequencies
    kernel = signal.windows.gaussian(kernel_width, std=(kernel_width) / 7.2)

    data_smooth = np.convolve(np.abs(channel_data - channel_data.mean()), kernel, 'same') / kernel.sum()

    # TODO: make threshold dependent on noise levels between events
    idxs_high = np.where(data_smooth > event_th)[0]  # indices where sound was ON
    
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
    periods = np.array(periods)  # these are all pairs of sample indices

    # convert samples indices into actual times in seconds, with the zero as the start of ephys recording
    # start of ephys recording is NOT the same as the start of ADC recording
    adc_timestamps = np.load(adc_ts_file)
    ephys_timestamps = np.load(ephys_ts_file)
    p_begs = adc_timestamps[periods[:, 0]] - ephys_timestamps[0]  # start of the period in seconds
    p_ends = adc_timestamps[periods[:, 1]] - ephys_timestamps[0]  # end of the period in seconds
    period_times = np.column_stack([p_begs, p_ends])

    events_exp = np.loadtxt(events_file, skiprows=1, delimiter=',')
    events_csv = np.loadtxt(sounds_file, skiprows=1, delimiter=',')
    events_csv[:, 0] = events_csv[:, 0] - events_exp[0][0]

    # first shift by the delay of the ephys start
    adc_to_ephys_shift = adc_timestamps[0] - ephys_timestamps[0]
    shift = events_csv[0][0] - period_times[0][0] + adc_to_ephys_shift
    events_csv[:, 0] = events_csv[:, 0] - shift

    # next linearly correct the drift by finding a time diff between
    # the last ADC pulse and the last logged one, assuming they are still the closest events
    t_last_ADC = period_times[-1][0]
    t_last_sev = events_csv[np.abs(t_last_ADC - events_csv[:, 0]).argmin()][0]
    drift = t_last_ADC - t_last_sev
    events_csv[:, 0] = events_csv[:, 0] + np.arange(len(events_csv)) * drift/len(events_csv)

    return period_times, events_csv


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