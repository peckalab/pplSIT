import numpy as np
import random


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