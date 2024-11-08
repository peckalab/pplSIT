from scipy import signal
from scipy.interpolate import interp1d
import numpy as np


def instantaneous_rate(spiketrain, time_bins, bin_size=10, k_width=70):
    """
    OLD legacy version! Use the one below!!

    spiketrain - an array of spike times in seconds
    time_bins  - an array of times to have instantaneous rate for
    """
    # add boundaries to match session length
    times = np.concatenate([np.array([0]), spiketrain, np.array([time_bins[-1]])])

    # histogram of spike times
    spikes_count, _ = np.histogram(times, bins=int(len(time_bins)/bin_size))

    # convolve with gaussian kernel for smoothing
    kernel = signal.gaussian(k_width, std=(k_width) / 7.2)
    i_rate = np.convolve(spikes_count/(1.0/bin_size), kernel, 'same') / kernel.sum()
    
    # interpolate to match experimental timeline
    lin = interp1d(np.linspace(0, time_bins[-1], len(i_rate)), i_rate)
    return lin(time_bins)


def inst_rate(spiketrain, time_bins, k_width=700):
    """
    spiketrain - an array of spike times in seconds
    time_bins  - an array of times to have instantaneous rate for
    """
    spikes_count, _ = np.histogram(spiketrain, bins=time_bins)

    # convolve with gaussian kernel for smoothing
    kernel = signal.gaussian(k_width, std=(k_width) / 7.2)
    i_rate = np.convolve(spikes_count/(1.0/np.diff(time_bins).mean()), kernel, 'same') / kernel.sum()

    return np.concatenate([i_rate, [i_rate[-1]]])  # double last value to match timeline


def spike_idxs(spiketrain, time_bins):
    """
    spiketrain - an array of spike times in seconds
    time_bins  - an array of times to have spike indices to
    """
    s_rate_tl = round(1.0 / np.diff(time_bins).mean())
    spiking_idxs = spiketrain * s_rate_tl
    spiking_idxs = spiking_idxs[spiking_idxs < len(time_bins)]
    return spiking_idxs.astype(np.int32)


def get_shuffled(spiketrain):
    # shuffle spike times preserving inter-spike intervals
    ISIs = np.diff(spiketrain)
    np.random.shuffle(ISIs)
    return np.concatenate([[spiketrain[0]], spiketrain[0] + np.cumsum(ISIs)])


def smooth_gaussian(data, k_width):
    kernel  = signal.gaussian(k_width, std=(k_width) / 7.2)
    return np.convolve(data, kernel, 'same') / kernel.sum()


def smooth_rectangular(data, width_in_bins):
    kernel = np.ones(width_in_bins)
    return np.convolve(data, kernel, 'same') / kernel.sum()