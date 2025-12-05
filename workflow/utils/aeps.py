import os
import numpy as np
from scipy import signal
from scipy import stats
from scipy.ndimage import convolve1d


AEP_metrics_lims = {
    'A1': {
        'P1': [10, 23],  # min - max to the left of min
        'N1': [23, 68],  # AUC
        'P2': [68, 100], # AUC
        'P3': [100, 200] # AUC
    },
    'PPC': {
        'P0': [5, 10], # min - max to the left of min
        'N0': [10, 15], # max - min to the left of max
        'P1': [15, 35], # min - max to the left of min
        'N1': [35, 70], # max - min to the left of max
        'P3': [70, 160] # AUC - consider 75 - 200!
    },
    'HPC': {
        'P0': [5, 10], # min - max to the left of min
        'N0': [10, 15], # max - min to the left of max
        'P1': [15, 35], # min - max to the left of min
        'N1': [35, 70], # max - min to the left of max
        'P3': [70, 160] # min - max to the left of min
    }
}

# 0 - min - max to the left of min, 1 - max - min to the left of max, 2 - AUC
AEP_metrics_methods = {
    'A1':  {'P1': 2, 'N1': 2, 'P2': 2, 'P3': 2},
    'PPC': {'P0': 0, 'N0': 1, 'P1': 0, 'N1': 1, 'P3': 2},
    'HPC': {'P0': 0, 'N0': 1, 'P1': 0, 'N1': 1, 'P3': 0},
}

# remove outliers (use > ~3SD?)
outlier_lims = {
    'A1': 5000,
    'PPC': 1500,
    'HPC': 1500
}

def compute_metric(aeps, method, t_l, t_r, k_width=20):
    metric = np.zeros(len(aeps))    
    if method == 0:  # min - max to the left of min
        min_idxs = aeps[:, t_l:t_r].argmin(axis=1)
        for i, idx in enumerate(min_idxs):
            if idx > 0:
                m_max = aeps[i][t_l:t_l + idx].max()
                m_min = aeps[i][t_l + idx]
                metric[i] = m_min - m_max
    elif method == 1:
        max_idxs = aeps[:, t_l:t_r].argmax(axis=1)
        for i, idx in enumerate(max_idxs):
            if idx > 0:
                m_min = aeps[i][t_l:t_l + idx].min()
                m_max = aeps[i][t_l + idx]
                metric[i] = m_max - m_min
    else:
        metric = aeps[:, t_l:t_r].sum(axis=1)
        
    # remove outliers?
    #m_mean, m_std = metric.mean(), metric.std()
    #metric[metric > m_mean + 3*m_std] = m_mean + 3*m_std
    #metric[metric < m_mean - 3*m_std] = m_mean - 3*m_std
    
    # z-scored
    metric_z = stats.zscore(metric)
    
    # smoothed
    kernel = signal.gaussian(k_width, std=(k_width) / 7.2)
    metric_smooth = np.convolve(metric_z, kernel, 'same') / kernel.sum()
    
    return metric, metric_smooth


def nanmean_filter1d(x, size, axis=0):
    # mask = 1 if finite else 0
    m = np.isfinite(x).astype(float)

    # replace NaN with 0
    x_filled = np.nan_to_num(x, nan=0.0)

    kernel = np.ones(size)

    # Sum of values
    x_sum = convolve1d(x_filled, kernel, axis=axis, mode='nearest')

    # Count of finite values
    m_sum = convolve1d(m, kernel, axis=axis, mode='nearest')

    # Avoid divide-by-zero
    out = x_sum / np.maximum(m_sum, 1e-8)
    out[m_sum == 0] = np.nan

    return out


def interp_nan_1d(x):
    x = np.asarray(x, float)
    nans = np.isnan(x)
    if not np.any(nans):
        return x

    t = np.arange(len(x))
    x[nans] = np.interp(t[nans], t[~nans], x[~nans])
    return x