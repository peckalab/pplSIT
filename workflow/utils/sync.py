import numpy as np
from scipy import signal
from utils.neurosuite import XMLHero, DatHero
import os


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
    adc_timestamps = np.load(adc_ts_file)
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
        if len(adc_timestamps) < 2:
            raise ValueError("Cannot infer ADC sample rate from fewer than 2 timestamps")
        dt = np.median(np.diff(adc_timestamps))
        if dt <= 0:
            raise ValueError(f"Invalid ADC timestamp spacing: median dt={dt}")
        s_rate = 1.0 / dt

    # --- load ADC data robustly ---
    raw = np.memmap(adc_file, dtype=dtype, mode="r")
    if raw.size % ch_no != 0:
        raise ValueError(
            f"ADC file size ({raw.size} values) is not divisible by channel count {ch_no} "
            f"for file {adc_file}"
        )

    n_frames = raw.size // ch_no

    if n_frames != len(adc_timestamps):
        raise ValueError(
            f"ADC frame count mismatch: dat implies {n_frames} frames, "
            f"timestamps has {len(adc_timestamps)} entries"
        )

    data = raw.reshape(n_frames, ch_no)
    channel_data = np.asarray(data[:, channel], dtype=float)

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