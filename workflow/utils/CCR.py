import numpy as np

def xcorr_session_with_shuffle(evoked, sustained, cond_indices, fs_hz=4, max_lag_s=20,
                               n_shuffle=200, shuffle_mode='cshift', seed=0):
    rng = np.random.default_rng(seed)
    T = len(evoked)
    L = int(round(max_lag_s * fs_hz))
    lags = None

    def build_event_train(n, idxs):
        e = np.zeros(n, float); e[np.asarray(idxs, int)] = 1.0; return e

    def lagged_pearson(x, y, L):
        n = len(x)
        lags = np.arange(-L, L+1)
        r = np.empty_like(lags, float)
        for i, lag in enumerate(lags):
            if lag >= 0:
                xs, ys = x[:n-lag], y[lag:]
            else:
                xs, ys = x[-lag:], y[:n+lag]
            m = np.isfinite(xs) & np.isfinite(ys)
            if m.sum() < 3: r[i] = np.nan; continue
            xs = (xs[m] - xs[m].mean()) / (xs[m].std(ddof=1)+1e-12)
            ys = (ys[m] - ys[m].mean()) / (ys[m].std(ddof=1)+1e-12)
            r[i] = np.mean(xs*ys)
        return lags, r

    def shuffle_event_train(e):
        if shuffle_mode == 'cshift':
            # shift by >= max_lag to avoid contaminating near-zero lags
            min_shift = L+1
            shift = rng.integers(min_shift, len(e)-min_shift) if len(e) > 2*min_shift else rng.integers(1, len(e))
            return np.roll(e, shift)
        elif shuffle_mode == 'permute':
            return rng.permutation(e)
        else:
            raise ValueError("shuffle_mode must be 'cshift' or 'permute'.")

    out = {}
    for cond, idxs in cond_indices.items():
        e_train = build_event_train(T, idxs)

        lags_samp, rE = lagged_pearson(evoked,   e_train, L)
        _,          rS = lagged_pearson(sustained, e_train, L)
        if lags is None: lags = lags_samp

        # Shuffled nulls
        rE_null = np.empty((n_shuffle, len(lags)))
        rS_null = np.empty((n_shuffle, len(lags)))
        for k in range(n_shuffle):
            e_shuf = shuffle_event_train(e_train)
            _, rE_null[k] = lagged_pearson(evoked,   e_shuf, L)
            _, rS_null[k] = lagged_pearson(sustained, e_shuf, L)

        out[cond] = {
            'lags_s': lags / fs_hz,
            'r_evoked': rE,   'r_sustained': rS,
            'rE_null_mean': np.nanmean(rE_null, axis=0),
            'rE_null_lo':   np.nanpercentile(rE_null, 2.5, axis=0),
            'rE_null_hi':   np.nanpercentile(rE_null, 97.5, axis=0),
            'rS_null_mean': np.nanmean(rS_null, axis=0),
            'rS_null_lo':   np.nanpercentile(rS_null, 2.5, axis=0),
            'rS_null_hi':   np.nanpercentile(rS_null, 97.5, axis=0),
        }
    return out