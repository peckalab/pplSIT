from scipy import signal
from scipy.interpolate import interp1d
import numpy as np


def instantaneous_rate(spiketrain, time_bins, bin_size=10, k_width=70):
    """
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


def spike_idxs(spiketrain, time_bins):
    """
    spiketrain - an array of spike times in seconds
    time_bins  - an array of times to have spike indices to
    """
    s_rate_tl = round(1.0 / np.diff(time_bins).mean())
    spiking_idxs = spiketrain * s_rate_tl
    spiking_idxs = spiking_idxs[spiking_idxs < len(time_bins)]
    return spiking_idxs.astype(np.int32)
