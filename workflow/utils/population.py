import os
import numpy as np
import h5py
from scipy import signal
from scipy import stats
from sklearn import decomposition


def unit_response_matrix(session_path, electrodes, times_to_event=[15, 28, 73, 100]):
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
