"""Data preparation and raw segmentations, independent of timing criteria."""
import math
import numpy as np

def select_events(sound_events, window_s, recording_bounds, time_ranges=None):
    """Return eligibility and exclusion flags for complete, nonoverlapping windows.

    Flags: 1 outside recording; 2 overlaps next logged sound event; 4 outside
    configured time ranges. Counts are kept per condition by the I/O layer.
    Equal exposure is required by the reference's unweighted pooled-spike fit.
    """
    events = np.asarray(sound_events, dtype=float)
    if (events.ndim != 2 or events.shape[1] < 2 or not np.all(np.isfinite(events[:, :2]))
            or np.any(np.diff(events[:, 0]) <= 0)
            or np.any(events[:, 1] != np.floor(events[:, 1]))):
        raise ValueError('sound_events must have finite time/ID columns and increasing times')
    bounds = np.asarray(recording_bounds, dtype=float)
    if bounds.shape != (2,) or not np.all(np.isfinite(bounds)) or bounds[0] >= bounds[1]:
        raise ValueError('Invalid recording bounds')
    if not math.isfinite(window_s) or window_s <= 0:
        raise ValueError('window_s must be finite and positive')
    starts = events[:, 0]
    stops = starts + window_s
    flags = np.zeros(len(events), dtype=np.uint8)
    # Tolerance is only for floating-point equality, never for acquisition jitter.
    eps = 1e-10
    flags[(starts < bounds[0] - eps) | (stops > bounds[1] + eps)] |= 1
    flags[:-1][stops[:-1] > starts[1:] + eps] |= 2
    if time_ranges is not None:
        ranges = np.asarray(time_ranges, dtype=float)
        if (ranges.ndim != 2 or ranges.shape[1] != 2 or not len(ranges)
                or not np.all(np.isfinite(ranges)) or np.any(ranges[:, 1] <= ranges[:, 0])
                or np.any(ranges[1:, 0] < ranges[:-1, 1])):
            raise ValueError('time_ranges_s must be sorted, disjoint [start, stop] intervals')
        in_range = np.zeros(len(events), dtype=bool)
        for start, stop in ranges:
            in_range |= (starts >= start - eps) & (stops <= stop + eps)
        flags[~in_range] |= 4
    return flags == 0, flags


def align_spikes(spike_times, event_times, window_s):
    """Pool spikes in [onset, onset + window); never include the next endpoint."""
    spikes = np.asarray(spike_times, dtype=float)
    times = np.asarray(event_times, dtype=float)
    if spikes.ndim != 1 or not np.all(np.isfinite(spikes)):
        raise ValueError('spike_times must be a finite vector')
    if times.ndim != 1 or not np.all(np.isfinite(times)) or np.any(np.diff(times) <= 0):
        raise ValueError('event_times must be a finite, strictly increasing vector')
    spikes = np.sort(spikes)
    lo = np.searchsorted(spikes, times, side='left')
    hi = np.searchsorted(spikes, times + window_s, side='left')
    chunks = [spikes[a:b] - t for t, a, b in zip(times, lo, hi) if b > a]
    return np.concatenate(chunks) if chunks else np.empty(0, dtype=float)



def empty_segmentation(status, n_events, n_spikes, diagnostic=''):
    return dict(status=status, diagnostic=diagnostic, n_events=int(n_events),
                n_spikes=int(n_spikes), fits=[])

def prepare_input(aligned_spikes, n_events, settings):
    """Shared validation/QC for CPU fitting and GPU batch preparation."""
    a = np.asarray(aligned_spikes, dtype=float)
    if (a.ndim != 1 or not np.all(np.isfinite(a))
            or np.any(a < 0) or np.any(a >= settings.window_s)):
        raise ValueError('Aligned spikes must lie in [0, window_s)')
    if isinstance(n_events, bool) or not isinstance(n_events, (int, np.integer)) or n_events < 0:
        raise ValueError('n_events must be a nonnegative integer')
    if n_events < settings.min_events:
        return a, empty_segmentation('insufficient_events', n_events, len(a))
    if len(a) < settings.min_spikes or len(np.unique(a)) < 2:
        return a, empty_segmentation('insufficient_spikes', n_events, len(a))
    if settings.max_spikes is not None and len(a) > settings.max_spikes:
        return a, empty_segmentation('spike_limit_exceeded', n_events, len(a))
    return a, None



def iter_segmentations(a, settings, backend=None):
    """Yield unfiltered raw fits; callers choose whether to consume every level."""
    if backend is None:
        from hepstats.modeling import bayesian_blocks
        backend = bayesian_blocks
    wrapped = np.concatenate((a - settings.window_s, a, a + settings.window_s))
    for p0 in settings.p0_values:
        try:
            with np.errstate(divide='raise', invalid='raise', over='raise'):
                raw = np.asarray(backend(wrapped, p0=p0), dtype=float)
        except FloatingPointError as exc:
            yield dict(p0=float(p0), status='backend_numerical_failure',
                       diagnostic=str(exc), raw_edges_s=np.empty(0))
            continue
        if (raw.ndim != 1 or len(raw) < 2 or not np.all(np.isfinite(raw))
                or np.any(np.diff(raw) <= 0)):
            raise ValueError('Bayesian-block backend returned invalid edges')
        yield dict(p0=float(p0), status='fitted', diagnostic='', raw_edges_s=raw)


def fit_spikes(aligned_spikes, n_events, settings, backend=None):
    """Cache every configured penalty, regardless of downstream stopping criteria."""
    a, unavailable = prepare_input(aligned_spikes, n_events, settings)
    if unavailable is not None:
        return unavailable
    result = empty_segmentation('fitted', n_events, len(a))
    result['fits'] = list(iter_segmentations(a, settings, backend))
    return result
