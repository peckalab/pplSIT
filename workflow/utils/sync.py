import numpy as np
from scipy import signal
from utils.neurosuite import XMLHero, DatHero


def get_sound_events_from_ephys(dat_file, xml_file, sounds_file, events_file, channel, event_th=200, ipi=0.25):
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
            events_synced[idx_ev][0] = pulse[0]

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