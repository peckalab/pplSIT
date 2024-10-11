import os
import numpy as np
import h5py
import json
from scipy import signal
from scipy import stats
from sklearn import decomposition
from utils.spiketrain import get_shuffled


def unit_activity_matrix(meta_file, units_file, electrodes, bin_size=0.01, shuffle=False):
    # electrodes - a list of electrodes to take into account, like [1, 2]
    with h5py.File(meta_file, 'r') as f:
        sound_events = np.array(f['processed']['sound_events'])
        cfg = json.loads(f['processed'].attrs['parameters'])

    spike_times = {}
    with h5py.File(units_file, 'r') as f:
        unit_names = [x for x in f if int(x.split('-')[0]) in electrodes]
    with h5py.File(units_file, 'r') as f:
        for unit_name in unit_names:
            spike_times[unit_name]  = np.array(f[unit_name]['spike_times'])

    # create binning first. Can't just do even bins
    # because sound events are not exactly spaced
    latency  = cfg['sound']['latency']  # seconds
    bins_per_event = int(latency / bin_size)
    event_bins = np.arange(0, bins_per_event) * bin_size

    # alternative way to solve odd sound latencies:
    # don't care about the last 5 bins because sometimes sound events are not equally spaced
    #event_bins = np.arange(0, bins_per_event - 5) * bin_size  # 10 ms bins
    #event_bins = np.concatenate([event_bins, event_bins[-1] + bin_size + np.arange(5) * 0.001])  # add 5 1 ms very little bins

    onsets = np.zeros([len(sound_events), bins_per_event])
    onsets[:, 0] = sound_events[:, 0]
    for i, bin_val in enumerate(event_bins):
        onsets[:, i] = sound_events[:, 0] + bin_val
        
    # correction for wrong delays - if normal binning used
    idxs_jitter = np.where(np.diff(onsets[:, 0]) < 0.241)[0]
    for idx in idxs_jitter:
        onsets[idx] = onsets[idx][0] + np.arange(bins_per_event)*0.0000001  # basically remove this pulse from analysis
    bins = onsets.flatten()  
    bins = np.concatenate([bins, [onsets[-1][0] + 0.25]])  # hack to have correct bins

    np.where(np.diff(bins) <= 0)[0]
    assert (np.diff(bins) > 0).all()  # just to make sure no overlaps

    # create activity matrix with "uneven" (almost even) bins
    unit_mx = np.zeros([len(sound_events) * bins_per_event, len(unit_names)])
    for i, unit_name in enumerate(unit_names):
        if shuffle:
            spiketrain = get_shuffled(spike_times[unit_name])
        else:
            spiketrain = spike_times[unit_name]
        unit_mx[:, i] = np.histogram(spiketrain, bins=bins)[0]
        #unit_mx[:, i] = stats.zscore(unit_mx[:, i])
        
    return bins, unit_mx.T  # units x time


def unit_response_matrix(session_path, electrodes, times_to_event=[15, 28, 73, 100]):
    """
    Legacy code for the unit activity matrix. Still can be used for e.g.
    specific periodic binning.
    """
    meta_file  = os.path.join(session_path, 'meta.h5')
    units_file = os.path.join(session_path, 'units.h5')

    with h5py.File(meta_file, 'r') as f:
        sound_events = np.array(f['processed']['sound_events'])

    spike_times = {}
    with h5py.File(units_file, 'r') as f:
        unit_names = [x for x in f if int(x.split('-')[0]) in electrodes]
    with h5py.File(units_file, 'r') as f:
        for unit_name in unit_names:
            spike_times[unit_name]  = np.array(f[unit_name]['spike_times'])

    # create binning first
    bins_per_event = len(times_to_event) + 1
    onsets = np.zeros([len(sound_events), bins_per_event])
    onsets[:, 0] = sound_events[:, 0]
    for i, m_lims in enumerate(times_to_event):
        onsets[:, i+1] = sound_events[:, 0] + m_lims/1000

    # correction for wrong delays
    idxs_jitter = np.where(np.diff(onsets[:, 0]) < 0.2)[0]
    for idx in idxs_jitter:
        onsets[idx] = onsets[idx][0] + np.arange(len(times_to_event) + 1)*0.000001  
    bins = onsets.flatten()  
    bins = np.concatenate([bins, [onsets[-1][0] + 0.25]])  # hack to have correct bins

    # create activity matrix with uneven bins
    unit_mx = np.zeros([len(sound_events) * bins_per_event, len(unit_names)])
    for i, unit_name in enumerate(unit_names):
        unit_mx[:, i] = np.histogram(spike_times[unit_name], bins=bins)[0]

    return bins, unit_mx


def activity_at_phase(s_path, phase=4, electrodes=[1, 2], do_pca=False, k_width=30):
    # by default = spontaneous activity, phase 4 (max 4)
    times_to_event = [15, 28, 73, 100]
    #times_to_event = [0, 50, 100, 150]

    bins, unit_mx = unit_response_matrix(s_path, electrodes, times_to_event)
    resp_at_phase = unit_mx[phase::len(times_to_event) + 1]
    unit_act_matrix = resp_at_phase.T
    for u, unit_data in enumerate(unit_act_matrix):
        unit_act_matrix[u] = stats.zscore(unit_data)
    resp_at_phase = unit_act_matrix.T

    if do_pca:
        pca = decomposition.PCA(n_components=3)
        pca.fit(resp_at_phase)
        X = pca.transform(resp_at_phase)
        pop_act = X[:, 0]  # PC1 score
    else:
        pop_act = resp_at_phase.mean(axis=1)  # or just a sum

    # smooth
    if k_width is not None:
        kernel  = signal.gaussian(k_width, std=(k_width) / 7.2)
        pop_act = np.convolve(pop_act, kernel, 'same') / kernel.sum()

    # filter slow oscillations
    sos = signal.butter(10, 0.001, fs=4, analog=False, btype='highpass', output='sos')
    pop_act = signal.sosfiltfilt(sos, pop_act)

    return pop_act


def pop_activity_phase_shifted(s_path, bins_bgr, bins_tgt, electrodes=[1, 2], do_pca=False, k_width=30):
    """
    Example binning (in ms) for diff durations:

    binning_mPFC_PCA = {
        50:  [8, 16, 65, 85],
        75:  [8, 16, 90, 110],
        100: [8, 16, 115, 135],
    }

    Returns: <sound_events_count> x <bin_count> matrix with the last column being pop activity
                between last bin and 250ms.
    """
    assert len(bins_bgr) == len(bins_tgt)

    # dependent on sound events!
    meta_file  = os.path.join(s_path, 'meta.h5')
    with h5py.File(meta_file, 'r') as f:
        events = np.array(f['processed']['sound_events'])
    idxs_bgr_ev = np.where(events[:, 1] == 1)[0]
    idxs_tgt_ev = np.where(events[:, 1] == 2)[0]

    bins, unit_mx_bgr = unit_response_matrix(s_path, electrodes, bins_bgr)
    bins, unit_mx_tgt = unit_response_matrix(s_path, electrodes, bins_tgt)

    w_mx = []
    for phase in range(len(bins_bgr)):
        bins_bgr_e = bins_bgr + [250]
        bins_tgt_e = bins_tgt + [250]

        # normalize by window size (important when different between BGR and TGT)
        resp_at_phase_bgr = unit_mx_bgr[phase + 1::len(bins_bgr) + 1] / ((bins_bgr_e[phase + 1] - bins_bgr_e[phase]) / 1000)
        resp_at_phase_tgt = unit_mx_tgt[phase + 1::len(bins_tgt) + 1] / ((bins_tgt_e[phase + 1] - bins_tgt_e[phase]) / 1000)

        # let no stimulus and noise episodes be binned as background
        resp_at_phase = resp_at_phase_bgr.copy()
        resp_at_phase[idxs_tgt_ev] = resp_at_phase_tgt[idxs_tgt_ev]

        unit_act_matrix = resp_at_phase.T
        for u, unit_data in enumerate(unit_act_matrix):
            unit_act_matrix[u] = stats.zscore(unit_data)
        resp_at_phase = unit_act_matrix.T

        if do_pca:
            pca = decomposition.PCA(n_components=3)
            pca.fit(resp_at_phase)
            X = pca.transform(resp_at_phase)
            pop_act = X[:, 0]  # PC1 score
        else:
            pop_act = resp_at_phase.mean(axis=1)  # or just a sum

        # smooth
        if k_width is not None:
            kernel  = signal.gaussian(k_width, std=(k_width) / 7.2)
            pop_act = np.convolve(pop_act, kernel, 'same') / kernel.sum()

        # filter slow oscillations
        sos = signal.butter(10, 0.001, fs=4, analog=False, btype='highpass', output='sos')
        pop_act = signal.sosfiltfilt(sos, pop_act)

        w_mx.append(stats.zscore(pop_act)) # stay in events space

    return np.column_stack(w_mx)