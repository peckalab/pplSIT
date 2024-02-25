"""
This module contains functions for computing multitaper spectrograms.

author: Arash Shahidi
email: A.Shahidi@campus.lmu.de
date: 29.01.2024

see also:
Elephant
https://elephant.readthedocs.io/en/latest/reference/_toctree/spectral/elephant.spectral.multitaper_psd.html

Prerau Lab
https://github.com/preraulab/multitaper_toolbox/blob/master/python/multitaper_spectrogram_python.py
"""

import numpy as np
import matplotlib.pyplot as plt
import pyfftw
from multiprocessing import cpu_count
from scipy import signal
from typing import Iterator, Tuple
from numpy.typing import ArrayLike

# %% utils
class SlidingWindow:
    """
    A class to create and manage sliding windows for time series data.

    This class facilitates the generation of sliding windows for a given dataset,
    allowing for the calculation of windowed metrics, such as in spectral analysis.
    It supports specifying window size, step size, and overlap ratio.

    Parameters
    ----------
    window_size : int
        The size of each window in terms of number of data points.
    window_step : int, optional
        The number of data points to step forward for the next window.
        If not provided, it is calculated based on the overlap ratio.
    overlap_ratio : float, optional
        The ratio of overlap between consecutive windows.
        If not provided, defaults to 0.5 (50% overlap).

    Attributes
    ----------
    window_size : int
        Size of the sliding window.
    window_step : int
        Step size between consecutive windows.

    Methods
    -------
    overlap()
        Returns the number of overlapping data points between consecutive windows.

    overlap_ratio()
        Returns the ratio of overlap between consecutive windows.

    overlap_ratio_to_window_step(overlap_ratio)
        Sets the window step size based on the provided overlap ratio.

    window_onsets(size)
        Computes the start indices of windows for a given dataset size.

    window_centers(size)
        Computes the center indices of windows for a given dataset size.

    window_centers_nwin(nwin)
        Computes the center indices for a specified number of windows.

    window_ends(size)
        Computes the end indices of windows for a given dataset size.

    window_indices(size)
        Generates the indices for each window for a given dataset size.

    window_view(X, axes=-1)
        Applies the sliding window to the input data along the specified axis.

    window_center_idx
        Returns the index of the center of the window.

    window_timesegs(size)
        Returns the start and end times of each window segment for a given size.

    window_generator(nchunks, size)
        Yields chunks of window indices for processing large datasets in parts.
    """
    def __init__(self, window_size: int, window_step: int = None, overlap_ratio: float = None):
        """
        Initialize the SlidingWindow object.

        Parameters
        ----------
        window_size : int
            The size of each window in terms of number of data points.
        window_step : int, optional
            The number of data points to step forward for the next window.
            If not provided, it is calculated based on the overlap ratio.
        overlap_ratio : float, optional
            The ratio of overlap between consecutive windows.
            If not provided, defaults to 0.5 (50% overlap).
        """
        self.window_size = window_size
        if window_step is not None:
            self.window_step = window_step
        else:
            self.overlap_ratio_to_window_step(0.5)
            
        if overlap_ratio is not None:
            self.overlap_ratio_to_window_step(overlap_ratio)

    @property
    def overlap(self) -> int:
        """
        Calculate the number of overlapping data points between consecutive windows.

        Returns
        -------
        int
            The number of overlapping data points.
        """
        return self.window_size - self.window_step

    @property
    def overlap_ratio(self) -> float:
        """
        Calculate the ratio of overlap between consecutive windows.

        Returns
        -------
        float
            The ratio of overlap as a decimal.
        """
        return self.overlap / self.window_size

    def overlap_ratio_to_window_step(self, overlap_ratio: float) -> None:
        """
        Set the window step size based on the provided overlap ratio.

        Parameters
        ----------
        overlap_ratio : float
            The desired ratio of overlap between windows.
        """
        self.window_step = self.window_size - int(overlap_ratio * self.window_size)

    def window_onsets(self, size: int) -> np.ndarray:
        """
        Compute the start indices of windows for a given dataset size.

        Parameters
        ----------
        size : int
            The size of the dataset.

        Returns
        -------
        np.ndarray
            An array of window start indices.
        """
        return np.arange(self.window_size-1, size, self.window_step)-(self.window_size-1)

    def window_centers(self, size: int) -> np.ndarray:
        """
        Compute the center indices of windows for a given dataset size.

        Parameters
        ----------
        size : int
            The size of the dataset.

        Returns
        -------
        np.ndarray
            An array of window center indices.
        """
        return self.window_onsets(size) + self.window_center_idx
    
    def window_centers_nwin(self, nwin: int) -> np.ndarray:
        """
        Compute the center indices for a specified number of windows.

        Parameters
        ----------
        nwin : int
            The number of windows.

        Returns
        -------
        np.ndarray
            An array of window center indices for the specified number of windows.
        """
        return np.arange(nwin) * self.window_step + self.window_center_idx

    def window_ends(self, size: int) -> np.ndarray:
        """
        Compute the end indices of windows for a given dataset size.

        Parameters
        ----------
        size : int
            The size of the dataset.

        Returns
        -------
        np.ndarray
            An array of window end indices.
        """
        return np.arange(self.window_size-1, size, self.window_step)

    def window_indices(self, size: int) -> np.ndarray:
        """
        Generate the indices for each window for a given dataset size.

        Parameters
        ----------
        size : int
            The size of the dataset.

        Returns
        -------
        np.ndarray
            A 2D array of indices for each window.
        """
        return (np.arange(self.window_size).reshape(1,-1) + self.window_onsets(size).reshape(-1,1)).astype(int)
    
    def window_view(self, X: np.ndarray, axes: int = -1) -> np.ndarray:
        """
        Apply the sliding window to the input data along the specified axis.

        Parameters
        ----------
        X : np.ndarray
            The input data array.
        axes : int, default -1
            The axis along which to apply the sliding window.

        Returns
        -------
        np.ndarray
            The data array segmented into windows.
            x_win: (len(x)-win_size+1, win_size)
        """
        axes %= len(X.shape)
        X = np.moveaxis(X, axes, -1)
        sliding_win = self.window_indices(X.shape[-1])
        X_win = X[..., sliding_win]
        X_win = unflush_axes(X_win, 2, axes)
        return X_win

    @property
    def window_center_idx(self) -> float:
        """
        Calculate the index of the center of the window.

        Returns
        -------
        float
            The index of the center of the window.
        """
        return (self.window_size-1)/2

    def window_timesegs(self, size: int) -> np.ndarray:
        """
        Return the start and end times of each window segment for a given size.

        Parameters
        ----------
        size : int
            The size of the dataset.

        Returns
        -------
        np.ndarray
            A 2D array with the start and end indices of each window.
        """
        return np.array([self.window_onsets(size), self.window_ends(size)]).T

    def window_generator(self, nchunks: int, size: int) -> Iterator[np.ndarray]:
        """
        Yield chunks of window indices for processing large datasets in parts.

        Parameters
        ----------
        nchunks : int
            The number of chunks to divide the windows into.
        size : int
            The size of the dataset.

        Yields
        ------
        Iterator[np.ndarray]
            Chunks of window indices.
        """
        # shape (nwindow, window_size)
        # split nwindow into smaller chunks to process one by one
        chunks = np.array_split(self.window_indices(size), nchunks)
        for chunk in chunks:
            yield chunk

def unflush_axes(X: np.ndarray, num_axes: int, dst_axis: int) -> np.ndarray:
    """
    Move the last `num_axes` axes of an array `X` to a new position specified by `dst_axis`.

    Args:
        X (np.ndarray): Input array.
        num_axes (int): Number of axes to move.
        dst_axis (int): Destination starting axis where the axes will be moved.

    Returns:
        np.ndarray: Array with the last `num_axes` axes moved to the new position.

    Example:
        X = np.random.rand(10, 20, 30, 40, 50)
        moved_array = unflush_axes(X, num_axes=3, dst_axis=1)
        print(moved_array.shape)  # (10, 30, 40, 50, 20)
    """
    last_axes = last_naxes(num_axes)
    dst_axes = tuple(dst_axis + np.arange(num_axes))
    X = np.moveaxis(X, last_axes, dst_axes)
    return X

def last_naxes(num_axes: int) -> Tuple[int, ...]:
    """
    Generate a tuple of the last `num_axes` axes in reverse order.

    Args:
        num_axes (int): The number of last axes to generate.

    Returns:
        Tuple[int, ...]: A tuple of the last `num_axes` axes in reverse order.

    Example:
        num_axes = 3
        print(last_naxes(num_axes))  # (-3, -2, -1)
    """
    shape_range = np.arange(num_axes)
    last_axes = tuple(-(shape_range + 1)[::-1])
    return last_axes

def pad_along_axis(arr, pad_width, axis=-2):
    """
    Pads the given array along a specific axis.
    
    :param arr: The array to pad.
    :param pad_width: A tuple (pad_before, pad_after) specifying the width of padding.
    :param axis: The axis along which to pad.
    :return: Padded array.
    """
    # Handle negative axis
    if axis < 0:
        axis += arr.ndim

    # Ensure the axis is valid
    if axis >= arr.ndim:
        raise ValueError("Invalid axis for padding.")

    # Create padding configuration
    padding_config = [(0, 0)] * arr.ndim
    padding_config[axis] = pad_width

    return np.pad(arr, padding_config)

#%% multitaper functions
def multitaper_psd(y: ArrayLike,
                   axis: int = -1,
                   NW: float = 2,
                   nfft: int = None,
                   K_max: int = None,
                   fs: float = 1,
                   detrend: bool = False,
                   n_fft_threads: int = cpu_count()) -> (np.ndarray, np.ndarray):
    """
    Compute the multitaper power spectral density (PSD) of a signal.

    Parameters
    ----------
    y : ArrayLike
        Input signal array.
    axis : int, default -1
        The axis along which the Fourier transform is applied.
    NW : float, default 2
        The time-halfbandwidth product.
    nfft : int, optional
        The number of points for the FFT computation. If None, defaults to the length of the signal.
    K_max : int, optional
        The maximum number of Slepian tapers. Defaults to 2*NW-1 if None.
    fs : float, default 1
        The sampling frequency of the input signal.
    detrend : bool, default False
        If True, detrend the signal before PSD computation.
    n_fft_threads : int, default cpu_count()
        The number of threads to use for FFT computation.

    Returns
    -------
    psd : np.ndarray
        The multitaper estimate of the signal's PSD.
    f : np.ndarray
        The frequency axis of the PSD.
    """
    y = np.swapaxes(y, axis, -1) # swap axes so the samples axis is the last dimension
    N = y.shape[-1] # number of samples
    if K_max is None:
        K_max = int(2*NW - 1)
    assert (K_max > 0 and K_max < N/2), print('increase resolution `len(y)` or decrease  `NW`')
    
    if detrend:
        y = signal.detrend(y, axis=-1)

    win = signal.windows.dpss(N, NW, Kmax=K_max)
    tapered_y = np.expand_dims(y, axis=-1) * win.T # (..., sample_dim, taper_dim)
    
    if nfft is None:
        nfft = N
    if nfft>N:
        tapered_y = pad_along_axis(tapered_y, (0, nfft - N), axis=-2)

    # psd, f = periodogram(tapered_y, fs=fs, freq_range=freq_range, axis=-2, detrend=detrend)
    f = np.fft.rfftfreq(nfft, d=1/fs)
    tapered_y = np.moveaxis(tapered_y, -2, 0)
    shape_origin = tapered_y.shape[1:]
    tapered_y = tapered_y.reshape(tapered_y.shape[0], -1)
    psd_shape = (len(f), tapered_y.shape[-1])
    # psd_shape = tapered_y.shape[:-2] + f.shape + (tapered_y.shape[-1],)
    psd_ = pyfftw.zeros_aligned(psd_shape, dtype='complex128')
    pyfftw.FFTW(tapered_y, psd_, axes=(0,), threads=n_fft_threads).execute()
    # f, psd = sciperiodogram(tapered_y, fs=fs, axis=-2, detrend=False, scaling='spectrum')
    psd_ = np.abs(psd_)**2/N
    psd_ = psd_.reshape(-1, *shape_origin)
    psd_ = np.mean(psd_, axis=-1) # average across tapers # (nfreq, nch, nwin, ntaper)
    psd_ = np.moveaxis(psd_, 0, -1) # nch, nwin, nfreq
    psd_ = np.swapaxes(psd_, axis, -1) # swap axes to original shape of input
    return psd_, f

def mtspecgram(X: ArrayLike,
               window_size: int,
               window_step: int,
               fs: float,
               NW: float,
               nfft: int = None,
               detrend: bool = True) -> (np.ndarray, np.ndarray, np.ndarray):
    """
    Compute a multitaper spectrogram of the input signal.

    Parameters
    ----------
    X : ArrayLike
        Input signal array. Can be either 1D (time,) or 2D (channels, time).
    window_size : int
        The size of each window for the spectrogram.
    window_step : int
        The step size between consecutive windows.
    fs : float
        The sampling frequency of the input signal.
    NW : float
        The time-halfbandwidth product.
    nfft : int, optional
        The number of points for the FFT computation. If None, defaults to the window size.
    detrend : bool, default True
        If True, detrend each window before PSD computation.

    Returns
    -------
    x_spec : np.ndarray (nch, nwin, nfreq), nch=1 if input is 1D
        The multitaper spectrogram.
    freqs : np.ndarray
        The frequency axis of the spectrogram.
    times : np.ndarray
        The time axis of the spectrogram.
    """
    if X.ndim == 1:
        X = X[np.newaxis, :]
    dur = X.shape[-1]
    slider = SlidingWindow(window_size, window_step)
    x_window = slider.window_view(X)
    times = slider.window_centers(dur) / fs
    x_spec, freqs = multitaper_psd(x_window, NW=NW, nfft=nfft, axis=-1, fs=fs, detrend=detrend)
    # ch, nwin, f -> ch, f, nwin
    x_spec = np.swapaxes(x_spec, -1, -2)
    return x_spec, freqs, times

# %% demo
import xarray as xr
def ar_simulate(freq, sample_rate, seconds, r=.95, noise_std=None, random_seed=None):
    """Create a simulated oscillation using an autoregressive filter.

    A simple filter is defined by direct pole placement and applied to white
    noise to generate a random signal with a defined oscillatory peak frequency
    that exhibits random variability frequency, amplitude and waveform.

    Parameters
    ----------
    freq : float
        Peak resonant frequency of the simulated filter.
    sample_rate : float
        Sampling frequency for the simulation
    seconds : float
        Number of seconds of data to simulate
    r : float (0 < r < 1)
        Pole magnitude of simulated autoregressive resonance.
    noise_std : float
        Scaling of optional noise to add to simulation. Scaling is relative to
        standard-deviation of the simulated data.
    random_seed : int
        Optional random seed generation

    Returns
    -------
    ndarray
        A simulated time course.

    Reference:
    [https://github.com/AJQuinn/emd-mirror/blob/master/emd/utils.py]
    """
    if random_seed is not None:
        np.random.seed(random_seed)

    if freq > 0:
        freq_rads = (2 * np.pi * freq) / sample_rate
        a1 = np.array([1, -2*r*np.cos(freq_rads), (r**2)])
    else:
        a1 = np.poly(r)

    num_samples = int(sample_rate * seconds)

    x = signal.filtfilt(1, a1, np.random.randn(1, num_samples)).T

    if noise_std is not None:
        noise = np.std(x)*noise_std*np.random.randn(1, num_samples).T
        x = x + noise

    if random_seed is not None:
        np.random.seed()  # restore defaults

    return x

if __name__=="__main__":
    # test multitaper_psd
    # ar_osillator 
    freq = 10
    fs = 1000
    T = 5
    y = ar_simulate(freq, fs, T, r=.95, noise_std=.1, random_seed=42)
    plt.plot(y)
    plt.show()

    mtspec, f = multitaper_psd(y.T[:, :1000], NW=2, fs=fs, nfft=1500, K_max=10, detrend=False)
    plt.plot(f, mtspec[0])
    plt.xlim([0, 20])
    plt.show()

    yspec, freq, time = mtspecgram(y.T, window_size=1024, window_step=500, fs=fs, NW=2, detrend=False)
    yspec = xr.DataArray(yspec, dims=['ch', 'freq', 'time'], coords={'ch': [0], 'freq': freq, 'time': time})
    yspec[0].sel(freq=slice(0,20)).plot.imshow()
    plt.show()

    yspec, freq, time = mtspecgram(y.T, window_size=1024, window_step=500, fs=fs, NW=2, nfft=1500, detrend=False)
    yspec = xr.DataArray(yspec, dims=['ch', 'freq', 'time'], coords={'ch': [0], 'freq': freq, 'time': time})
    yspec[0].sel(freq=slice(0,20)).plot.imshow()
    plt.show()

