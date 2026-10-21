"""Eight-bin, presentation-normalized tuning following code_07062021.

The bootstrap resamples whole sound presentations, including zero-spike ones.
Aggregating identical (bin, spike-count) pairs gives exactly the same sampling
distribution as resampling individual event indices, without a large event array.
"""
import hashlib
import json
from pathlib import Path
import warnings

import numpy as np

CENTERS = np.arange(8)*np.pi/4
PHASES = ('onset', 'offset', 'late')
TABLE = json.loads(Path(__file__).with_name('low_spike_thresholds.json').read_text())


def settings(value):
    result = dict(timing_condition='pooled', conditions=None, references=None,
                  unit_names=None, comparisons=None, n_bootstrap=1000, seed=0)
    unknown = set(value)-set(result)
    if unknown:
        raise ValueError(f'Unknown egocentric_tuning settings: {sorted(unknown)}')
    result.update(value)
    for key in ('n_bootstrap', 'seed'):
        v = result[key]
        if isinstance(v, bool) or not isinstance(v, int) or v < (2 if key == 'n_bootstrap' else 0):
            raise ValueError(f'{key} must be an integer >= {2 if key == "n_bootstrap" else 0}')
    for key in ('conditions', 'references', 'unit_names'):
        v = result[key]
        if v is not None and (not isinstance(v, list) or not v
                              or any(not isinstance(n, str) or not n or '/' in n for n in v)
                              or len(set(v)) != len(v)):
            raise ValueError(f'{key} must be null or a nonempty list of unique names')
    v = result['timing_condition']
    if v is not None and (not isinstance(v, str) or not v or '/' in v):
        raise ValueError('timing_condition must be null or a condition name')
    pairs = result['comparisons']
    if pairs is not None:
        if not isinstance(pairs, list) or any(
                not isinstance(p, list) or len(p) != 2 or p[0] == p[1]
                or any(not isinstance(n, str) or not n or '/' in n for n in p) for p in pairs):
            raise ValueError('comparisons must be null or a list of distinct-name pairs [A, B]')
        if len({tuple(p) for p in pairs}) != len(pairs):
            raise ValueError('Duplicate comparisons')
    return result


def phase_windows(onset, offset, duration, observation):
    """Intervals in seconds relative to sound onset; late includes pre-onset."""
    missing = {p: (np.empty((0, 2)), 'onset_unavailable') for p in PHASES}
    if not np.isfinite(onset) or onset < 0:
        return missing
    if not np.isfinite(duration) or duration <= 0:
        return {p: (np.empty((0, 2)), 'duration_unavailable') for p in PHASES}
    detected = np.isfinite(offset)
    end = offset if detected else onset+duration
    candidates = dict(onset=[[onset, end]], offset=[[end, end+duration]],
                      late=[[0., onset], [end+duration, observation]])
    result = {}
    for phase, intervals in candidates.items():
        a = np.asarray(intervals, dtype=float)
        if (not np.isfinite(observation) or observation <= 0
                or end <= onset or np.any(a[:, 1] < a[:, 0])
                or np.any(a < 0) or np.any(a > observation)):
            result[phase] = (np.empty((0, 2)), 'window_outside_observation')
        else:
            a = a[a[:, 1] > a[:, 0]]
            result[phase] = (a, ('detected_offset' if detected else 'duration_fallback')
                             if len(a) else 'zero_duration')
    return result


def count_spikes(spikes, events, intervals):
    """Half-open intervals [start, stop), consistent across adjacent phases."""
    counts = np.zeros(len(events), dtype=np.int64)
    for start, stop in intervals:
        counts += (np.searchsorted(spikes, events+stop, side='left')
                   - np.searchsorted(spikes, events+start, side='left'))
    return counts


def angle_bins(angles):
    if not np.all(np.isfinite(angles)):
        raise ValueError('Angular binning requires finite bearings')
    return np.floor(((np.rad2deg(angles)+22.5) % 360)/45).astype(int)


def vector_metrics(curves):
    curves = np.asarray(curves, dtype=float)
    values = np.nan_to_num(curves, nan=0.)
    total = values.sum(axis=-1, keepdims=True)
    norm = np.divide(values, total, out=np.full_like(values, np.nan), where=total > 0)
    directional = np.sum(norm*np.exp(1j*CENTERS), axis=-1)
    folded = np.sum(norm*np.exp(2j*CENTERS), axis=-1)
    return np.stack([np.abs(directional), np.abs(folded),
                     np.angle(directional) % (2*np.pi),
                     (np.angle(folded) % (2*np.pi))/2], axis=-1)


def stable_seed(seed, *names):
    encoded = json.dumps([seed, *names], separators=(',', ':')).encode()
    return int.from_bytes(hashlib.sha256(encoded).digest()[:8], 'little')


def bootstrap_curves(bins, counts, repetitions, seed):
    if len(bins) == 0:
        return np.full((repetitions, 8), np.nan)
    pairs, frequency = np.unique(np.column_stack([bins, counts]), axis=0, return_counts=True)
    draws = np.random.default_rng(seed).multinomial(len(bins), frequency/frequency.sum(),
                                                 size=repetitions)
    curves = np.zeros((repetitions, 8))
    for b in range(8):
        selected = pairs[:, 0] == b
        denominator = draws[:, selected].sum(axis=1)
        numerator = (draws[:, selected]*pairs[selected, 1]).sum(axis=1)
        # Reference bootstrap zero-fills unsampled bins before metric calculation.
        curves[:, b] = np.divide(numerator, denominator, out=np.zeros(repetitions),
                                where=denominator > 0)
    return curves


def percentile(values, levels, axis=0):
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', message='All-NaN slice encountered')
        return np.nanpercentile(values, levels, axis=axis)


def criteria(spikes, presentations, metrics, lower):
    """Return reference classification, preserving strict low-count lookup."""
    point_threshold = np.nan
    if spikes < 75:
        row = [r for r in TABLE['rows'] if r['No_spikes'] < spikes]
        thresholds = np.array([row[-1]['VS'], row[-1]['VS2fold']]) if row else np.ones(2)
        passed = np.asarray(lower) > thresholds
    else:
        point_threshold = .2
        thresholds = np.array([.15, .15])
        passed = (np.asarray(metrics[:2]) >= .2) & (np.asarray(lower) >= thresholds)
    if np.min(presentations) < 10:
        status = 'insufficient_angular_coverage'
    elif np.max(presentations)/np.min(presentations) > 10:
        status = 'asymmetric_angular_coverage'
    elif spikes <= 20:
        status = 'insufficient_spikes'
    else:
        status = 'eligible'
    selected = int(np.argmax(metrics[:2]))
    classification = ('directional', 'orientation')[selected] if (
        status == 'eligible' and passed[selected]) else ('none' if status == 'eligible' else 'unavailable')
    return dict(status=status, classification=classification,
                directional_pass=bool(passed[0] and status == 'eligible'),
                orientation_pass=bool(passed[1] and status == 'eligible'),
                directional_lower_threshold=float(thresholds[0]),
                orientation_lower_threshold=float(thresholds[1]),
                point_threshold=point_threshold)


def analyze(angles, counts, duration, repetitions, seed):
    counts = np.asarray(counts)
    if (len(angles) != len(counts) or np.any(counts < 0)
            or np.any(counts != np.floor(counts)) or not np.isfinite(duration) or duration <= 0):
        raise ValueError('Require matched nonnegative integer counts and positive phase duration')
    bins = angle_bins(angles)
    occupancy = np.bincount(bins, minlength=8)
    spike_hist = np.bincount(bins, weights=counts, minlength=8)
    curve = np.divide(spike_hist, occupancy, out=np.full(8, np.nan), where=occupancy > 0)
    metrics = vector_metrics(curve)
    boot = bootstrap_curves(bins, counts, repetitions, seed)
    boot_metrics = vector_metrics(boot)
    interval = percentile(boot_metrics[:, :2], [15.9, 50, 84.1])
    criterion = criteria(int(counts.sum()), occupancy, metrics, interval[0])
    axis = metrics[3]
    maximum = CENTERS[np.argmax(np.nan_to_num(curve, nan=0.))]
    delta = lambda a: abs((a-maximum+np.pi) % (2*np.pi)-np.pi)
    oriented = axis if delta(axis) < delta(axis+np.pi) else (axis+np.pi) % (2*np.pi)
    # Zero-strength vectors do not have meaningful preferred angles.
    direction = metrics[2] if metrics[0] > 1e-12 else np.nan
    axis = axis if metrics[1] > 1e-12 else np.nan
    oriented = oriented if metrics[1] > 1e-12 else np.nan
    ci = percentile(boot, [2.5, 97.5])
    ci[:, occupancy == 0] = np.nan
    result = dict(presentations=occupancy, spikes=spike_hist.astype(np.int64),
                  spikes_per_presentation=curve, rate_hz=curve/duration,
                  rate_ci95_hz=ci/duration, vector_strength_interval=interval,
                  n_events=len(counts), n_spikes=int(counts.sum()), duration_s=duration,
                  vector_strength=float(metrics[0]), orientation_strength=float(metrics[1]),
                  preferred_direction_rad=float(direction), preferred_axis_rad=float(axis),
                  preferred_orientation_rad=float(oriented),
                  bootstrap_valid_metrics=int(np.isfinite(boot_metrics[:, 0]).sum()),
                  bootstrap_seed=str(seed), **criterion)
    return result, boot/duration


def compare(a, b, boot_a, boot_b):
    """B minus A, independent-presentation bootstrap; no multiplicity claim."""
    delta = b['rate_hz']-a['rate_hz']
    boot_delta = boot_b-boot_a
    ci = percentile(boot_delta, [2.5, 97.5])
    missing = (a['presentations'] == 0) | (b['presentations'] == 0)
    ci[:, missing] = np.nan
    full = not np.any(missing)
    result = dict(rate_difference_hz=delta, rate_difference_ci95_hz=ci,
                  equal_angle_mean_difference_hz=float(np.mean(delta)) if full else np.nan,
                  equal_angle_mean_difference_ci95_hz=(
                      percentile(boot_delta.mean(axis=1), [2.5, 97.5]) if full else np.full(2, np.nan)),
                  vector_strength_difference=b['vector_strength']-a['vector_strength'],
                  orientation_strength_difference=b['orientation_strength']-a['orientation_strength'],
                  direction_shift_rad=(b['preferred_direction_rad']-a['preferred_direction_rad']
                                       + np.pi) % (2*np.pi)-np.pi,
                  axis_shift_rad=(b['preferred_axis_rad']-a['preferred_axis_rad']
                                  + np.pi/2) % np.pi-np.pi/2,
                  status='eligible' if a['status'] == b['status'] == 'eligible' else 'descriptive_only')
    return result
