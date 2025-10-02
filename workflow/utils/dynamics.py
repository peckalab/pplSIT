# Re-run after kernel reset
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import welch
from scipy.stats import entropy

def compute_temporal_dynamics(lfp_ptp, states, fs=4, window_sec=2, shuffle_n=100, state_names=None):
    """
    Analyze temporal dynamics of LFP_ptp per behavioral state using ACF, PSD, and sliding variance/entropy.

    Parameters:
    - lfp_ptp: 1D array of instantaneous metric (e.g., evoked LFP PTP), shape (n_pulses,)
    - states: 1D array of behavioral states, shape (n_pulses,)
    - fs: sampling frequency (in Hz) of the pulses (default 4 Hz)
    - window_sec: sliding window size in seconds
    - shuffle_n: number of shuffles for null model
    - state_names: optional list of strings to label each state
    """

    def autocorr(x):
        x = x - np.nanmean(x)
        result = np.correlate(x, x, mode='full')
        mid = len(result) // 2
        result = result[mid:] / result[mid]  # normalize
        return result

    def sliding_metrics(x, win_size):
        var_list, ent_list = [], []
        for i in range(len(x) - win_size + 1):
            window = x[i:i + win_size]
            if np.any(np.isnan(window)):
                var_list.append(np.nan)
                ent_list.append(np.nan)
                continue
            hist, _ = np.histogram(window, bins=10, density=True)
            ent_list.append(entropy(hist + 1e-10))  # avoid log(0)
            var_list.append(np.var(window))
        return np.array(var_list), np.array(ent_list)

    unique_states = np.unique(states)
    max_lag = 20  # max lag in time points for ACF
    win_size = int(window_sec * fs)

    fig, axes = plt.subplots(len(unique_states), 4, figsize=(20, 4 * len(unique_states)))
    if state_names is None:
        state_names = [f"State {s}" for s in unique_states]

    for i, s in enumerate(unique_states):
        idx = states == s
        signal = lfp_ptp[idx]
        time = np.arange(len(signal)) / fs

        # ACF
        acf_real = autocorr(signal)[:max_lag]
        acf_nulls = np.array([autocorr(np.random.permutation(signal))[:max_lag] for _ in range(shuffle_n)])
        acf_null_mean = np.mean(acf_nulls, axis=0)
        acf_null_std = np.std(acf_nulls, axis=0)

        # PSD
        f_psd, psd_real = welch(signal, fs=fs, nperseg=min(len(signal), 64))
        psd_nulls = np.array([welch(np.random.permutation(signal), fs=fs, nperseg=min(len(signal), 64))[1]
                              for _ in range(shuffle_n)])
        psd_null_mean = np.mean(psd_nulls, axis=0)
        psd_null_std = np.std(psd_nulls, axis=0)

        # Sliding window variance & entropy
        var_vals, ent_vals = sliding_metrics(signal, win_size)

        # === Plot ACF
        ax = axes[i, 0]
        lags = np.arange(max_lag) / fs
        ax.plot(lags, acf_real, label='Real', color='black')
        ax.fill_between(lags, acf_null_mean - acf_null_std, acf_null_mean + acf_null_std,
                        color='gray', alpha=0.3, label='Shuffled ±1STD')
        ax.set_title(f"{state_names[i]} - ACF")
        ax.set_xlabel("Lag (s)")
        ax.set_ylabel("ACF")
        ax.legend()

        # === Plot PSD
        ax = axes[i, 1]
        ax.semilogy(f_psd, psd_real, label='Real', color='black')
        ax.fill_between(f_psd, psd_null_mean - psd_null_std, psd_null_mean + psd_null_std,
                        color='gray', alpha=0.3, label='Shuffled ±1STD')
        ax.set_title(f"{state_names[i]} - PSD")
        ax.set_xlabel("Frequency (Hz)")
        ax.set_ylabel("Power")
        ax.legend()

        # === Plot Sliding Variance
        ax = axes[i, 2]
        ax.plot(time[:len(var_vals)], var_vals, color='royalblue')
        ax.set_title(f"{state_names[i]} - Sliding Variance")
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Variance")

        # === Plot Sliding Entropy
        ax = axes[i, 3]
        ax.plot(time[:len(ent_vals)], ent_vals, color='firebrick')
        ax.set_title(f"{state_names[i]} - Sliding Entropy")
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Entropy")

    plt.tight_layout()
    plt.show()
