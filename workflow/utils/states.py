import numpy as np
import h5py
import json
import os

from scipy import signal
from utils.behavior import get_idxs_as_periods


def get_state_as_periods(s_path, sound_state, loc_state, att_state, min_dur, speed_th=0.04, strip_l=2, strip_r=1):
    """
    s_path - session path

    Options:
    - sound_state: BGR, TGT, SIL
    - loc_state: RUN, STA or None
    - att state: AL, PH or None

    min_dur - minimum duration in sound pulses.
    """
    meta_file = os.path.join(s_path, 'meta.h5')

    with h5py.File(meta_file, 'r') as f:
        tl = np.array(f['processed']['timeline'])
        sound_events = np.array(f['processed']['sound_events'])
        cfg = json.loads(f['processed'].attrs['parameters'])
        tgt_mx = np.array(f['processed']['target_matrix'])

    speed_ev = tl[sound_events[:, 2].astype(np.int32)][:, 3]

    idxs_dict = {
        'BGR': np.where(sound_events[:, 1] == 1)[0],
        'TGT': np.where(sound_events[:, 1] == 2)[0],
        'SIL': np.where(sound_events[:, 1] == 0)[0],
        'DI1': np.where(sound_events[:, 1] == 3)[0],
        'DI2': np.where(sound_events[:, 1] == 4)[0],
        'RUN': np.where(speed_ev > speed_th)[0],
        'STA': np.where(speed_ev < speed_th)[0],
    }

    if att_state is not None:
        ensembles_file = os.path.join(s_path, 'analysis', 'ensembles.h5')
        if not os.path.exists(ensembles_file):
            raise FileNotFoundError('Ensembles file not found')
        with h5py.File(ensembles_file, 'r') as f:
            ens_ev = np.array(f['AL'])[sound_events[:, 2].astype(np.int32)]

        # apply optional smoothing, aimed to increase periods a bit
        width = 40  # it's actually 10 seconds
        kernel = signal.windows.gaussian(width, std=(width) / 7.2)
        ens_ev_smooth = np.convolve(ens_ev, kernel, 'same') / kernel.sum()

        idxs_dict['AL'] = np.where(ens_ev_smooth > 0)[0]
        idxs_dict['PH'] = np.where(ens_ev_smooth < 0)[0]

    # times of first tgt success pulses - just in case
    idxs_tgt_first_ev = tgt_mx[tgt_mx[:, 4] == 1][:, 0]
    tgt_first_t = sound_events[idxs_tgt_first_ev][:, 0]

    # success stays - just in case
    idxs_tgt_succ = []
    for tgt_rec in tgt_mx[tgt_mx[:, 4] == 1]:
        idxs_tgt_succ += list(np.arange(tgt_rec[0], tgt_rec[1] + 1))
    idxs_tgt_succ = np.array(idxs_tgt_succ)

    # get indices
    idxs_selected = idxs_dict[sound_state]

    if loc_state is not None:
        idxs_selected = np.intersect1d(idxs_selected, idxs_dict[loc_state])

    if att_state is not None:
        idxs_selected = np.intersect1d(idxs_selected, idxs_dict[att_state])

    # convert to periods for selected condition
    periods_loc_ev = get_idxs_as_periods(idxs_selected)
    periods_filt = periods_loc_ev[np.where(np.diff(periods_loc_ev, axis=1) > min_dur)[0]]

    # don't consider first X and last Y pulses as boundary condition
    periods_filt_ex = []
    for per in periods_filt:
        per_ex = [per[0] + strip_l, per[1] - strip_r]
        if per_ex[1] - per_ex[0] > 0:  # at least one pulse
            periods_filt_ex.append(per_ex)
    periods_filt_ex = np.array(periods_filt_ex)

    long_pers_xy = np.zeros([len(periods_filt_ex), 4], dtype=object)
    for k, ls_rec in enumerate(periods_filt_ex):
        idx_tl_s = int(sound_events[ls_rec[0]][2])
        idx_tl_e = int(sound_events[ls_rec[1]][2])
        x_pos = tl[np.arange(idx_tl_s, idx_tl_e)][:, 1]
        y_pos = tl[np.arange(idx_tl_s, idx_tl_e)][:, 2]

        long_pers_xy[k] = [int(ls_rec[0]), int(ls_rec[1]), float(x_pos.mean()), float(y_pos.mean())]

    # back to indices
    idxs_long_pers_ev = []
    for ls in long_pers_xy:
        idxs_long_pers_ev += list(np.arange(int(ls[0]), int(ls[1])))
    idxs_long_pers_ev = np.array(idxs_long_pers_ev) # sound_events[idxs_long_pers_ev][:, 0]

    # first is a periods count x 4 matrix, sampled at sound events. idx start, idx end, X, Y
    # second is the same but just a list of pulse indices
    return long_pers_xy, idxs_long_pers_ev