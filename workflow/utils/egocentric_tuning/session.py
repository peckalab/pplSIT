"""Join geometry, reusable response timing, and spikes by original event index."""
import csv
import json
from pathlib import Path
import time

import h5py
import numpy as np

from utils.egocentric_geometry.session import file_record
from .core import (CENTERS, PHASES, TABLE, analyze, compare, count_spikes,
                   phase_windows, settings, stable_seed)


def selection(requested, available, name):
    result = list(available) if requested is None else list(requested)
    missing = set(result)-set(available)
    if not result or missing:
        raise ValueError(f'Unavailable {name}: {sorted(missing)}; available: {list(available)}')
    return result


def write_result(group, result):
    for key, value in result.items():
        if isinstance(value, np.ndarray):
            group.create_dataset(key, data=value)
        else:
            group.attrs[key] = value


def scalar_row(result):
    return {k: v for k, v in result.items() if not isinstance(v, np.ndarray)}


def write_csv(path, rows, required):
    fields = list(dict.fromkeys([*required, *(k for row in rows for k in row)]))
    with Path(path).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fields)
        writer.writeheader()
        writer.writerows(rows)


def check_join(geometry, responses, units_path):
    if (geometry.attrs.get('artifact_type') != 'egocentric_sound_geometry'
            or not geometry.attrs.get('complete', False)):
        raise ValueError('Require a complete egocentric geometry artifact')
    if responses.attrs.get('method') != 'amaro_07062021_bayesian_blocks':
        raise ValueError('Require Bayesian-block response timing from code_07062021')
    events = responses['sound_events'][:]
    if not np.array_equal(events, geometry['source_sound_events'][:], equal_nan=True):
        raise ValueError('Geometry and response sound events differ; regenerate with consistent synchronization')
    indices = geometry['event_indices'][:]
    if (indices.dtype.kind not in 'iu' or len(indices) == 0
            or np.any(np.diff(indices) <= 0) or indices[0] < 0 or indices[-1] >= len(events)
            or not np.array_equal(geometry['event_times_s'][:], events[indices, 0])
            or not np.array_equal(geometry['sound_ids'][:], events[indices, 1])):
        raise ValueError('Geometry event-index mapping is inconsistent')
    records = json.loads(responses.attrs['input_files'])
    record = [r for r in records if Path(r['path']).name == 'units.h5']
    stat = Path(units_path).stat()
    if (len(record) != 1 or record[0]['size_bytes'] != stat.st_size
            or record[0]['mtime_ns'] != stat.st_mtime_ns):
        raise ValueError('Units export differs from response timing provenance; regenerate responses '
                         'using these units (or use the original units export)')
    lookup = np.full(len(events), -1, dtype=np.int64)
    lookup[indices] = np.arange(len(indices))
    return events, lookup


def run_session(inputs, output_h5, output_csv, output_comparisons, config, session):
    started = time.monotonic()
    config = settings(config)
    outputs = [Path(p).resolve() for p in (output_h5, output_csv, output_comparisons)]
    if len(set(outputs)) != 3 or set(outputs) & {Path(p).resolve() for p in inputs.values()}:
        raise ValueError('Tuning outputs must be distinct from inputs and each other')
    rows, comparison_rows = [], []
    with (h5py.File(inputs['geometry'], 'r') as geometry,
          h5py.File(inputs['responses'], 'r') as responses,
          h5py.File(inputs['units'], 'r') as units):
        events, lookup = check_join(geometry, responses, inputs['units'])
        if geometry.attrs['session'] != session:
            raise ValueError('Geometry session differs from requested session')
        conditions = selection(config['conditions'], [
            c for c in json.loads(responses.attrs['condition_names']) if c != 'pooled'], 'conditions')
        references = selection(config['references'], json.loads(geometry.attrs['references']), 'references')
        names = selection(config['unit_names'], json.loads(responses.attrs['unit_names']), 'units')
        timing_condition = config['timing_condition']
        if timing_condition is not None and timing_condition not in responses['conditions']:
            raise ValueError(f'Missing timing condition {timing_condition}; enable include_pooled or select another')
        pairs = config['comparisons']
        if pairs is None:
            pairs = [conditions] if len(conditions) == 2 else []
        for pair in pairs:
            selection(pair, conditions, 'comparison conditions')
        window = float(json.loads(responses.attrs['effective_settings'])['window_s'])
        condition_data = {}
        for condition in conditions:
            group = responses['conditions'][condition]
            indices = group['event_indices'][:]
            if (indices.dtype.kind not in 'iu' or np.any(np.diff(indices) <= 0)
                    or np.any(indices < 0) or np.any(indices >= len(events))
                    or not np.array_equal(group['event_times_s'][:], events[indices, 0])):
                raise ValueError(f'Inconsistent response event mapping: {condition}')
            mapped = lookup[indices]
            if np.any(mapped < 0):
                raise ValueError(f'{condition} includes events without positive-source geometry rows')
            data = dict(indices=indices, times=events[indices, 0], mapped=mapped)
            data['references'] = {}
            for ref in references:
                g = geometry['events'][ref]
                angles = g['bearing_rad'][:][mapped]
                valid = (g['geometry_flags'][:][mapped] == 0) & np.isfinite(angles)
                data['references'][ref] = (angles, valid)
            condition_data[condition] = data
        for a, b in pairs:
            if np.intersect1d(condition_data[a]['indices'], condition_data[b]['indices']).size:
                raise ValueError(f'Comparison {a}, {b} shares events; independent-condition bootstrap requires disjoint events')
        for path in outputs:
            path.parent.mkdir(parents=True, exist_ok=True)
        source_files = [Path(__file__), Path(__file__).with_name('core.py'),
                        Path(__file__).with_name('low_spike_thresholds.json')]
        with h5py.File(output_h5, 'w') as out:
            out.attrs.update(
                artifact_type='egocentric_neural_tuning', schema_version=1, complete=False,
                method='amaro_07062021_centered_presentation_bootstrap', session=session,
                config=json.dumps(config), conditions=json.dumps(conditions), references=json.dumps(references),
                unit_names=json.dumps(names), comparisons=json.dumps(pairs),
                input_files=json.dumps({k: file_record(p) for k, p in inputs.items()}),
                implementation_sha256=json.dumps({p.name: file_record(p)['sha256'] for p in source_files}),
                reference=json.dumps({k: v for k, v in TABLE.items() if k != 'rows'}),
                stimulus_catalog=responses.attrs['stimulus_catalog'],
                observation_window_s=window,
                interval_convention='[start, stop), seconds from sound onset',
                angle_convention='radians; 0=in front; positive=counterclockwise in arena x/y',
                event_selection='response-eligible events intersected with mathematically valid geometry',
                tracking_quality_policy='unchanged upstream synchronized tracking; no added quality filters',
                comparison_convention='B minus A; independent presentation bootstrap; pointwise 95% CIs, uncorrected',
                temporal_adaptation='actual sound duration; late includes pre-onset and post-offset intervals',
                bootstrap_method='multinomial counts of identical (angle-bin, spike-count) presentation pairs',
                classification='Amaro criteria, not an auditory responsiveness test or a p-value',
                package_versions=json.dumps(dict(numpy=np.__version__, h5py=h5py.__version__)))
            out.create_dataset('angular_centers_rad', data=CENTERS)
            for condition, data in condition_data.items():
                g = out.create_group(f'conditions/{condition}')
                original = responses['conditions'][condition]
                for key in ('label', 'sound_ids', 'stimuli', 'duration_s', 'duration_status'):
                    if key in original.attrs:
                        g.attrs[key] = original.attrs[key]
                g.attrs['n_response_eligible_events'] = len(data['indices'])
                g.create_dataset('event_indices', data=data['indices'])
                g.create_dataset('event_times_s', data=data['times'])
                for ref, (angles, valid) in data['references'].items():
                    rg = g.create_group('geometry/'+ref)
                    rg.create_dataset('bearing_rad', data=angles)
                    rg.create_dataset('valid', data=valid)
            for unit in names:
                spikes = units[unit]['spike_times'][:]
                if spikes.ndim != 1 or not np.all(np.isfinite(spikes)) or np.any(np.diff(spikes) < 0):
                    raise ValueError(f'Unit {unit} spikes must be finite and sorted')
                computed = {}
                for condition, data in condition_data.items():
                    source = responses['conditions'][condition]
                    original = source['units'][unit].attrs
                    whole_counts = count_spikes(spikes, data['times'], [[0., window]])
                    if int(whole_counts.sum()) != original['n_spikes']:
                        raise ValueError(f'{unit}/{condition}: spikes disagree with stored response counts')
                    tc = timing_condition or condition
                    timing = responses['conditions'][tc]['units'][unit].attrs
                    onset = float(timing.get('onset_latency_s', np.nan))
                    offset = float(timing.get('offset_latency_s', np.nan))
                    if timing.get('onset_status') != 'detected':
                        onset = np.nan
                    if timing.get('offset_status') != 'detected':
                        offset = np.nan
                    windows = phase_windows(onset, offset, float(source.attrs['duration_s']), window)
                    ug = out.create_group(f'conditions/{condition}/units/{unit}')
                    ug.attrs.update(timing_condition=tc, onset_latency_s=onset, offset_latency_s=offset,
                                    onset_status=timing['onset_status'], offset_status=timing['offset_status'])
                    for phase, (intervals, timing_status) in windows.items():
                        pg = ug.create_group(phase)
                        pg.create_dataset('intervals_s', data=intervals)
                        pg.attrs['timing_status'] = timing_status
                        duration = float(np.sum(intervals[:, 1]-intervals[:, 0]))
                        counts = count_spikes(spikes, data['times'], intervals) if duration > 0 else None
                        if counts is not None:
                            pg.create_dataset('event_spike_counts', data=counts, compression='gzip')
                        for ref, (angles, valid) in data['references'].items():
                            rg = pg.create_group(ref)
                            info = dict(unit=unit, condition=condition, reference=ref, phase=phase,
                                        timing_condition=tc, timing_status=timing_status,
                                        n_events_excluded_geometry=int((~valid).sum()))
                            if counts is None:
                                result = dict(status='timing_unavailable', classification='unavailable',
                                              duration_s=duration, n_events=int(valid.sum()))
                            else:
                                result, boot = analyze(angles[valid], counts[valid], duration,
                                                       config['n_bootstrap'], stable_seed(
                                                           config['seed'], unit, condition, ref, phase))
                                computed[condition, ref, phase] = (result, boot)
                            write_result(rg, result)
                            rows.append({**info, **scalar_row(result)})
                for pair_index, (a, b) in enumerate(pairs):
                    for ref in references:
                        for phase in PHASES:
                            if (a, ref, phase) in computed and (b, ref, phase) in computed:
                                va, ba = computed[a, ref, phase]
                                vb, bb = computed[b, ref, phase]
                                result = compare(va, vb, ba, bb)
                            else:
                                result = dict(status='timing_unavailable')
                            g = out.create_group(f'comparisons/{pair_index}/{unit}/{ref}/{phase}')
                            write_result(g, result)
                            comparison_rows.append(dict(unit=unit, reference=ref, phase=phase,
                                                        condition_a=a, condition_b=b, **scalar_row(result)))
                print(f'Tuning {unit}: {len(computed)} curves', flush=True)
            write_csv(output_csv, rows, ['unit', 'condition', 'reference', 'phase', 'status'])
            write_csv(output_comparisons, comparison_rows,
                      ['unit', 'reference', 'phase', 'condition_a', 'condition_b', 'status'])
            out.attrs['elapsed_seconds'] = time.monotonic()-started
            out.attrs['complete'] = True
