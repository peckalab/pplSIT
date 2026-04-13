import numpy as np


def get_event_periods(tl, event_type):
    # event_type: -1, 0, 1, 2 (noise, silence, background, target)
    # returns: periods in seconds of event_type

    idxs_events = np.where(tl[:, 6] == event_type)[0]

    # no events of this type
    if len(idxs_events) == 0:
        return np.zeros((0, 2), dtype=float)

    # find breaks between contiguous stretches
    idxs_to_idxs = np.where(np.diff(idxs_events) > 1)[0]

    # one single contiguous period
    if len(idxs_to_idxs) == 0:
        return np.array([[tl[idxs_events[0], 0], tl[idxs_events[-1], 0]]], dtype=float)

    # multiple periods
    periods = np.zeros((len(idxs_to_idxs) + 1, 2), dtype=np.int32)

    # first block
    periods[0] = np.array([0, idxs_to_idxs[0]])

    # middle blocks
    if len(idxs_to_idxs) > 1:
        periods[1:-1] = np.column_stack([idxs_to_idxs[:-1] + 1, idxs_to_idxs[1:]])

    # last block
    periods[-1] = np.array([idxs_to_idxs[-1] + 1, len(idxs_events) - 1])

    # convert from positions within idxs_events -> indices in tl -> times
    tl_idxs_mx = np.column_stack([idxs_events[periods[:, 0]], idxs_events[periods[:, 1]]])
    return np.column_stack([tl[tl_idxs_mx[:, 0], 0], tl[tl_idxs_mx[:, 1], 0]])


def get_sound_event_periods(sound_events, event_type):
    t_periods = []
    curr_period = []
    # if sound_events[0][1] == event_type:  # event starts with the first pulse
    #     curr_period.append(sound_events[0][0])
    for i in range(len(sound_events) - 1):  # always starts with BGR, so ignore first pulse
        if sound_events[i-1][1] != event_type and sound_events[i][1] == event_type:  # start of the period
            curr_period.append(sound_events[i][0])
        if sound_events[i+1][1] != event_type and sound_events[i][1] == event_type:  # end of the period
            # !!! time of the FIRST PULSE AFTER period end
            curr_period.append(sound_events[i+1][0])

            t_periods.append(curr_period)
            curr_period = []
    return t_periods


def get_sound_event_period_idxs(sound_events, event_type):
    # TODO: make DRY, union with the function above
    t_periods = []
    curr_period = []
    for i in range(len(sound_events) - 1):  # always starts with BGR, so ignore first pulse
        if sound_events[i-1][1] != event_type and sound_events[i][1] == event_type:  # start of the period
            curr_period.append(i)
        if sound_events[i+1][1] != event_type and sound_events[i][1] == event_type:  # end of the period
            # !!! time of the FIRST PULSE AFTER period end
            curr_period.append(i+1)

            t_periods.append(curr_period)
            curr_period = []
    return t_periods