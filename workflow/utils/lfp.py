import numpy as np
import random

from scipy.signal import windows
from scipy.signal import butter, filtfilt, hilbert, iirnotch


def detect_artifacts_amp(x, fs, win_ms=50, step_ms=10, z_thresh=8):
    """
    x : 1D LFP signal
    fs : sampling rate (Hz)
    win_ms : window length in ms
    step_ms: step between windows in ms
    z_thresh : threshold in robust z units
    """
    win = int(win_ms * fs / 1000)
    step = int(step_ms * fs / 1000)

    # compute RMS in sliding windows
    rms_vals = []
    idx = []
    for start in range(0, len(x) - win, step):
        seg = x[start:start+win]
        rms_vals.append(np.sqrt(np.mean(seg**2)))
        idx.append(start + win//2)
    rms_vals = np.array(rms_vals)
    idx = np.array(idx)

    # robust baseline: median & MAD
    med = np.median(rms_vals)
    mad = np.median(np.abs(rms_vals - med)) + 1e-12  # avoid /0

    z = np.abs(rms_vals - med) / mad  # "robust z-score"

    # windows that are extreme
    bad = z > z_thresh

    # convert window centers → sample mask
    mask = np.zeros_like(x, dtype=bool)
    for c, b in zip(idx, bad):
        if b:
            start = max(0, c - win//2)
            end   = min(len(x), c + win//2)
            mask[start:end] = True

    return mask, z, idx


def detect_artifacts_hf(x, fs, hp_hz=150, win_ms=50, z_thresh=8):
    # high-pass
    b, a = butter(4, hp_hz/(fs/2), btype='high')
    x_hf = filtfilt(b, a, x)

    # envelope
    env = np.abs(hilbert(x_hf))

    # smooth envelope over window
    win = int(win_ms * fs / 1000)
    kernel = np.ones(win) / win
    env_smooth = np.convolve(env, kernel, mode='same')

    # robust z on envelope
    med = np.median(env_smooth)
    mad = np.median(np.abs(env_smooth - med)) + 1e-12
    z = (env_smooth - med) / mad

    mask = z > z_thresh
    return mask, z


def clean_mask(mask, fs, min_dur_ms=50, merge_gap_ms=50):
    """
    mask: 1D bool array, True = artifact
    fs:   sampling rate (Hz)
    min_dur_ms: minimum duration of an artifact region to keep
    merge_gap_ms: merge two artifact regions if the gap between them
                  is <= this duration
    """
    mask = np.asarray(mask, bool)
    n = mask.size

    # convert ms → samples
    min_len   = int(round(min_dur_ms * fs / 1000.0))
    merge_gap = int(round(merge_gap_ms * fs / 1000.0))

    idx = np.flatnonzero(mask)
    if idx.size == 0:
        return np.zeros_like(mask, bool)

    # ---- 1) Find contiguous runs in the original mask ----
    # idx are sorted; breaks where gap > 1 sample
    breaks = np.where(np.diff(idx) > 1)[0]
    starts = np.r_[idx[0], idx[breaks + 1]]
    ends   = np.r_[idx[breaks], idx[-1]]

    # ---- 2) Merge runs whose gaps are <= merge_gap ----
    merged_starts = []
    merged_ends   = []

    cur_start = starts[0]
    cur_end   = ends[0]

    for s, e in zip(starts[1:], ends[1:]):
        # gap between current merged block and next block
        gap = s - cur_end - 1
        if gap <= merge_gap:
            # extend current block
            cur_end = e
        else:
            # close current block
            merged_starts.append(cur_start)
            merged_ends.append(cur_end)
            # start new
            cur_start, cur_end = s, e

    # flush last block
    merged_starts.append(cur_start)
    merged_ends.append(cur_end)

    # ---- 3) Apply minimum duration and build cleaned mask ----
    cleaned = np.zeros_like(mask, bool)
    for s, e in zip(merged_starts, merged_ends):
        if (e - s + 1) >= max(min_len, 1):
            cleaned[s:e+1] = True

    return cleaned


def global_rms_trace(lfp):
    """
    lfp: (n_samples, n_channels)
    returns: 1D array, RMS across channels for each time point
    """
    return np.sqrt(np.mean(lfp**2, axis=1))


def compute_bootstrapped_baseline(stationary_segments, running_segments, n_bootstraps=100, seed=42, ci_percentile=95):
    """
    Computes a bootstrapped balanced baseline across stationary and running segments.

    Parameters:
    - stationary_segments, running_segments: list of LFP segments (each segment is time x channels)
    - n_bootstraps: number of bootstrap repetitions
    - seed: random seed for reproducibility
    - ci_percentile: confidence interval width (default: 95 for 2.5% to 97.5%)

    Returns:
    - baseline_mean: mean baseline per channel (averaged across bootstraps)
    - baseline_std: std per channel across bootstraps
    - baseline_ci_lower: lower bound of confidence interval per channel
    - baseline_ci_upper: upper bound of confidence interval per channel
    """
    #random.seed(seed)
    #np.random.seed(seed)

    min_len = min(len(stationary_segments), len(running_segments))
    if min_len == 0:
        raise ValueError("Not enough segments in one or both groups for bootstrapping.")

    bootstrapped_means = []

    for _ in range(n_bootstraps):
        sampled_stat = random.sample(stationary_segments, min_len)
        sampled_run = random.sample(running_segments, min_len)
        all_segments = sampled_stat + sampled_run

        # Compute mean LFP of each segment (averaged over time)
        segment_means = [seg.mean(axis=0) for seg in all_segments]
        bootstrapped_means.append(np.mean(segment_means, axis=0))

    bootstrapped_means = np.stack(bootstrapped_means, axis=0)  # shape: (n_bootstraps, channels)

    baseline_mean = np.mean(bootstrapped_means, axis=0)
    baseline_std  = np.std(bootstrapped_means, axis=0)

    # Confidence intervals
    lower_percentile = (100 - ci_percentile) / 2
    upper_percentile = 100 - lower_percentile

    baseline_ci_lower = np.percentile(bootstrapped_means, lower_percentile, axis=0)
    baseline_ci_upper = np.percentile(bootstrapped_means, upper_percentile, axis=0)

    return baseline_mean, baseline_std, baseline_ci_lower, baseline_ci_upper


def notch_50hz(lfp, fs, freq=50.0, q=30.0):
    """
    lfp: (n_samples, n_channels) or (n_samples,)
    fs: sampling rate (Hz)
    freq: line frequency
    q: quality factor (higher = narrower notch)
    """
    b, a = iirnotch(freq, Q=q, fs=fs)
    # axis=0 assumes time is axis 0; adjust if needed
    lfp_filt = filtfilt(b, a, lfp, axis=0)
    return lfp_filt

def lowpass_lfp(lfp, fs, cutoff=40.0, order=4):
    nyq = fs / 2.0
    b, a = butter(order, cutoff/nyq, btype='low')
    return filtfilt(b, a, lfp, axis=0)

def highpass_lfp(lfp, fs, cutoff=40.0, order=4):
    nyq = fs / 2.0
    b, a = butter(order, cutoff/nyq, btype='high')
    return filtfilt(b, a, lfp, axis=0)