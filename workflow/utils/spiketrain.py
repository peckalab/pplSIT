from scipy import signal
from scipy.stats import pearsonr
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


def unit_response_metrics_per_condition(spike_times, stimuli, conditions, response_window=(0, 0.25), \
                                        bin_size=0.01, evoked_bins=(0, 12), baseline_bins=(15, 25)):
    """
    Compute extended response metrics for one unit (spiketrain) per each condition given (stimuli).

    Parameters:
    - spike_times: 1D array of spike times (in seconds)
    - stimuli: stimulus times (in seconds)
    - conditions: dict of condition indices to stimuli array, like {'idxs_bgr_sta': [1, 2, 3, ..., 9485], }
    - response_window: tuple, response window in seconds (e.g., 0–0.25)
    - bin_size: bin size for PSTH in seconds
    - total_window: total window for trial-aligned spiking (used for PSTH and Fano)

    Returns:
    - Dictionary: condition_id -> metrics
    """
    total_window = response_window[1] - response_window[0]
    n_bins = int(total_window / bin_size)
    bin_edges = np.linspace(0, total_window, n_bins + 1)
    metrics = {}
    #conditions = np.unique(stimuli[:, 1])

    for cond_id, cond_idxs in conditions.items():
        if len(cond_idxs) < 5:
            continue

        #stim_times = stimuli[stimuli[:, 1] == cond][:, 0]
        stim_times = stimuli[cond_idxs]
        n_trials = len(stim_times)

        aligned_spikes = []
        spike_counts = []
        first_spike_latencies = []

        # Create trial-by-trial binned spike trains
        trial_psths = np.zeros((n_trials, n_bins))
        for i, stim in enumerate(stim_times):
            trial_spikes = spike_times[(spike_times >= stim) & (spike_times < stim + total_window)] - stim
            aligned_spikes.append(trial_spikes)
            trial_hist, _ = np.histogram(trial_spikes, bins=bin_edges)
            trial_psths[i, :] = trial_hist / bin_size  # convert to rate
            spike_counts.append(len(trial_spikes[(trial_spikes >= response_window[0]) & (trial_spikes <= response_window[1])]))

            # latency
            spikes_in_window = trial_spikes[(trial_spikes >= response_window[0]) & (trial_spikes <= response_window[1])]
            if len(spikes_in_window) > 0:
                first_spike_latencies.append(spikes_in_window[0])
        
        # Metrics
        mean_rate = np.mean(spike_counts) / (response_window[1] - response_window[0])
        psth_avg = np.mean(trial_psths, axis=0)
        ptp = np.ptp(psth_avg)
        ptp_evoked   = np.ptp(psth_avg[evoked_bins[0]:evoked_bins[1]])
        ptp_baseline = np.ptp(psth_avg[baseline_bins[0]:baseline_bins[1]])
        std_baseline = np.std(psth_avg[baseline_bins[0]:baseline_bins[1]])
        rms = np.sqrt(np.mean(psth_avg**2))
        rms_norm = rms / (mean_rate + 1e-6)
        latency_ms = np.mean(first_spike_latencies) * 1000 if len(first_spike_latencies) > 0 else np.nan
        normalized_ptp = ptp / (mean_rate + 1e-3)
        snr_base_ptp = ptp / (ptp_baseline + 1e-3)
        snr_base_std = ptp / (std_baseline + 1e-3)

        # Response reliability (split-half corr)
        split1 = np.mean(trial_psths[:n_trials//2], axis=0)
        split2 = np.mean(trial_psths[n_trials//2:], axis=0)
        if np.std(split1) > 0 and np.std(split2) > 0:
            reliability = pearsonr(split1, split2)[0]
        else:
            reliability = np.nan

        # Temporal jitter (std of first spike latencies)
        jitter = np.std(first_spike_latencies) * 1000 if len(first_spike_latencies) > 1 else np.nan

        # Fano Factor
        fano = np.var(spike_counts) / (np.mean(spike_counts) + 1e-6)

        if ptp_baseline < 1e-6:
            ptp_ratio = np.inf if ptp_evoked > 0 else 0
        else:
            ptp_ratio = ptp_evoked / ptp_baseline

        metrics[cond_id] = {
            'mean_rate': mean_rate,
            'ptp': ptp,
            'rms': rms,
            'rms_norm': rms_norm,
            'latency_ms': latency_ms,
            'reliability': reliability,
            'jitter_ms': jitter,
            'fano_factor': fano,
            #'psth': psth_avg,
            'normalized_ptp': normalized_ptp,
            'ptp_evoked': ptp_evoked,
            'ptp_baseline': ptp_baseline,
            'ptp_ratio': ptp_ratio,
            'snr_base_ptp': snr_base_ptp,
            'snr_base_std': snr_base_std
        }

    return metrics


def compute_sttc_per_condition(spike_trains_dict, t_window=0.02, total_duration=None, event_conditions=None, conditions_to_use=None):
    """
    Compute STTC separately for each condition.

    Parameters:
    - spike_trains_dict: dict {unit_id: spike_times in seconds}
    - t_window: time window for tiling (default 20ms)
    - total_duration: total duration of the entire recording (optional)
    - event_conditions: 2D array (n_events x 2): [event_time, condition_id]
    - conditions_to_use: list or set of condition_ids to compute (optional)

    Returns:
    - sttc_results: dict {condition_id: {(unit1, unit2): STTC}}
    """

    def proportion_spikes_near(spikes_A, spikes_B, delta):
        if len(spikes_A) == 0 or len(spikes_B) == 0:
            return 0.0
        spikes_A = np.sort(spikes_A)
        spikes_B = np.sort(spikes_B)
        count = 0
        idx_B = 0
        for spike in spikes_A:
            while idx_B < len(spikes_B) and spikes_B[idx_B] < spike - delta:
                idx_B += 1
            if idx_B < len(spikes_B) and abs(spikes_B[idx_B] - spike) <= delta:
                count += 1
        return count / len(spikes_A)

    def compute_T(spikes, delta, duration):
        if len(spikes) == 0:
            return 0.0
        intervals = []
        for s in spikes:
            start = max(0, s - delta)
            end = min(duration, s + delta)
            intervals.append((start, end))
        intervals.sort()
        merged = []
        for start, end in intervals:
            if not merged or start > merged[-1][1]:
                merged.append([start, end])
            else:
                merged[-1][1] = max(merged[-1][1], end)
        total_time = sum(e - s for s, e in merged)
        return total_time / duration

    # === Determine conditions ===
    if event_conditions is None:
        raise ValueError("event_conditions must be provided to compute STTC per condition.")
    unique_conditions = np.unique(event_conditions[:, 1])
    if conditions_to_use is not None:
        unique_conditions = [c for c in unique_conditions if c in conditions_to_use]

    sttc_results = {}
    for cond in unique_conditions:
        # Determine time limits for this condition
        times_cond = event_conditions[event_conditions[:, 1] == cond][:, 0]
        if len(times_cond) == 0:
            continue
        cond_start = times_cond.min()
        cond_end = times_cond.max() + 0.25  # assume 250ms response window
        cond_duration = cond_end - cond_start if total_duration is None else total_duration

        # Slice spikes for this condition
        spike_trains_cond = {
            unit: spikes[(spikes >= cond_start) & (spikes <= cond_end)] - cond_start
            for unit, spikes in spike_trains_dict.items()
        }

        # Compute STTC matrix for this condition
        units = list(spike_trains_cond.keys())
        sttc_cond = {}
        for unit1, unit2 in combinations(units, 2):
            spikes_A = spike_trains_cond[unit1]
            spikes_B = spike_trains_cond[unit2]

            PA = proportion_spikes_near(spikes_A, spikes_B, t_window)
            PB = proportion_spikes_near(spikes_B, spikes_A, t_window)

            TA = compute_T(spikes_A, t_window, cond_duration)
            TB = compute_T(spikes_B, t_window, cond_duration)

            denom_A = 1 - PA * TB
            denom_B = 1 - PB * TA
            sttc = 0.5 * ((PA - TB) / denom_A + (PB - TA) / denom_B) if denom_A > 0 and denom_B > 0 else np.nan
            sttc_cond[(unit1, unit2)] = sttc

        sttc_results[cond] = sttc_cond

    return sttc_results