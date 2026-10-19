"""Amaro boundary selection, independent of the numerical segmentation backend."""
import math
import numpy as np

from .settings import Settings
from .segmentation import align_spikes, select_events, prepare_input, iter_segmentations

from .metadata import METHOD, METHOD_VERSION

def _empty_result(status, n_events, n_spikes):
    return dict(status=status, diagnostic='', is_modulated=-1, n_events=int(n_events), n_spikes=int(n_spikes),
                onset_status='not_estimated', onset_latency_s=np.nan, modulation_p0=np.nan,
                onset_edge_p0=np.nan, offset_status='not_estimated', offset_latency_s=np.nan,
                offset_p0=np.nan, short_latency_adjusted=False, fits=[])



def estimate_timing(aligned_spikes, n_events, duration_s, settings, backend=None):
    """Convenience API preserving lazy reference fitting for a single case."""
    a, unavailable = prepare_input(aligned_spikes, n_events, settings)
    if unavailable is not None:
        return select_timing([], n_events, len(a), duration_s, settings,
                             status=unavailable['status'], diagnostic=unavailable['diagnostic'])
    return select_timing(iter_segmentations(a, settings, backend), n_events, len(a),
                         duration_s, settings)


def select_timing(raw_fits, n_events, n_spikes, duration_s, settings,
                  status='fitted', diagnostic=''):
    """Select response boundaries without fitting or requiring aligned spikes.

    p0 is segmentation regularization, not a responsiveness-test p-value.
    Increases and decreases can both supply boundaries; no baseline test is used.
    """
    if status != 'fitted':
        if status not in ('insufficient_events', 'insufficient_spikes',
                          'spike_limit_exceeded', 'backend_numerical_failure'):
            raise ValueError('Unknown segmentation status')
        result = _empty_result(status, n_events, n_spikes)
        result['diagnostic'] = diagnostic
        return result
    result = _empty_result('analyzed', n_events, n_spikes)
    result.update(is_modulated=0, onset_status='not_modulated', offset_status='no_onset')
    previous = None
    latency = np.nan
    last_changes = None
    iterator = iter(raw_fits)
    for p0 in settings.p0_values:
        try:
            saved = next(iterator)
        except StopIteration as exc:
            raise ValueError('Incomplete segmentation cache: missing penalty level') from exc
        if saved['p0'] != p0:
            raise ValueError('Segmentation penalty levels do not match settings')
        if saved['status'] == 'backend_numerical_failure':
            unavailable = _empty_result('backend_numerical_failure', n_events, n_spikes)
            unavailable['diagnostic'] = saved['diagnostic']
            return unavailable
        if saved['status'] != 'fitted':
            raise ValueError('Unknown segmentation status')
        raw = np.asarray(saved['raw_edges_s'], dtype=float)
        if (raw.ndim != 1 or len(raw) < 2 or not np.all(np.isfinite(raw))
                or np.any(np.diff(raw) <= 0)):
            raise ValueError('Invalid cached Bayesian-block edges')
        # Dated reference: remove the RIGHT edge of each <1 ms block, once.
        edges = np.delete(raw, np.flatnonzero(np.diff(raw) < settings.min_block_s) + 1)
        fit = dict(p0=float(p0), raw_edges_s=raw, edges_s=edges)
        result['fits'].append(fit)
        changes = np.flatnonzero((edges > settings.edge_margin_s)
                                & (edges < settings.window_s - settings.edge_margin_s))
        if not len(changes):
            break
        result.update(is_modulated=1, modulation_p0=float(p0))
        # Preserve the reference's previous-successful-level latency selection.
        if previous is None:
            latency = edges[changes[0]]
            result['onset_edge_p0'] = float(p0)
        else:
            latency = previous['edges_s'][last_changes[0]]
            result['onset_edge_p0'] = previous['p0']
        previous, last_changes = fit, changes

    if result['is_modulated'] == 0:
        return result
    if latency < settings.short_latency_s:
        result['short_latency_adjusted'] = True
        first = last_changes[0]
        candidate = first + 1 if latency == previous['edges_s'][first] else first
        latency = previous['edges_s'][candidate] if candidate < len(previous['edges_s']) else np.nan
        result['onset_edge_p0'] = previous['p0']
    if not (settings.edge_margin_s < latency < settings.window_s - settings.edge_margin_s):
        result.update(onset_status='reference_boundary_outside_window', offset_status='no_valid_onset')
        return result
    result.update(onset_latency_s=float(latency), onset_status='detected')
    if duration_s is None or not np.isfinite(duration_s) or duration_s <= 0:
        result['offset_status'] = 'duration_unavailable'
        return result
    low, high = latency + duration_s + np.asarray(settings.offset_search_relative_to_duration_s)
    if low <= latency or high >= settings.window_s:
        result['offset_status'] = 'search_outside_window'
        return result
    result['offset_status'] = 'not_detected'
    # Offset fitting reuses the exact fits already computed for onset selection.
    for fit in result['fits']:
        if fit['p0'] < result['modulation_p0']:
            break
        candidates = fit['edges_s'][(fit['edges_s'] > low) & (fit['edges_s'] < high)]
        if not len(candidates):
            break
        result.update(offset_status='detected', offset_latency_s=float(candidates[0]), offset_p0=fit['p0'])
    return result


def diagnostic_psth(aligned_spikes, n_events, settings):
    edges = np.linspace(0, settings.window_s, math.ceil(settings.window_s / settings.psth_bin_s) + 1)
    counts = np.histogram(aligned_spikes, edges)[0]
    rate = counts / (n_events * np.diff(edges)) if n_events else np.full(len(counts), np.nan)
    return edges, counts, rate
