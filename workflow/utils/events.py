import numpy as np


def get_event_periods(tl, event_type):
    # event_type: -1, 0, 1, 2 (noise, silence, background, target)
    # returns: periods in seconds of event_type
    idxs_events  = np.where(tl[:, 6] == event_type)[0]
    idxs_to_idxs = np.where(np.diff(idxs_events) > 1)[0]

    # periods - indices to TL where was silent
    periods       = np.zeros([len(idxs_to_idxs) + 1, 2])
    periods[0]    = np.array([0, idxs_to_idxs[0]])
    periods[1:-1] = np.column_stack([idxs_to_idxs[:-1] + 1, idxs_to_idxs[1:]])
    periods[-1]   = np.array([idxs_to_idxs[-1], len(idxs_events) - 1])
    periods       = periods.astype(np.int32)

    # convert to TL indices and then to times
    tl_idxs_mx = np.column_stack([idxs_events[periods[:, 0]], idxs_events[periods[:, 1]]])
    return np.column_stack([tl[tl_idxs_mx[:, 0]][:, 0], tl[tl_idxs_mx[:, 1]][:, 0]]) 


def get_sound_event_periods(sound_events, event_type):
    t_periods = []
    curr_period = []
    for i in range(len(sound_events) - 1):  # always starts with BGR, so ignore first pulse
        if sound_events[i-1][1] != event_type and sound_events[i][1] == event_type:  # start of the period
            curr_period.append(sound_events[i][0])
        if sound_events[i+1][1] != event_type and sound_events[i][1] == event_type:  # end of the period
            # !!! time of the FIRST PULSE AFTER period end
            curr_period.append(sound_events[i+1][0])

            t_periods.append(curr_period)
            curr_period = []
    return t_periods