"""Cheap response-boundary selection from a complete segmentation cache."""
import csv
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import time

import h5py
import numpy as np

from .settings import Settings, FIT_SETTING_KEYS, fitting_config, timing_config
from .metadata import CACHE_ARTIFACT, CACHE_SCHEMA_VERSION
from .bayesian_blocks import _empty_result, select_timing, diagnostic_psth
# Preserve the existing public helper import for callers.
from .fits_session import make_conditions, _input_record


def select_session(fits_path, output_h5, output_csv, config):
    """Select timing and rebin PSTHs without importing or invoking a fitting backend."""
    paths = [Path(p).resolve() for p in (fits_path, output_h5, output_csv)]
    if len(set(paths)) != 3:
        raise ValueError('Cache, response HDF5 and summary CSV paths must be different')
    started = time.perf_counter()
    with open(fits_path, 'rb') as stream:
        cache_hash = hashlib.file_digest(stream, 'sha256').hexdigest()
    summary = []
    with h5py.File(fits_path, 'r') as source:
        if (source.attrs.get('artifact_type') != CACHE_ARTIFACT
                or source.attrs.get('schema_version') != CACHE_SCHEMA_VERSION
                or not source.attrs.get('complete', False)):
            raise ValueError('A complete fits.h5 segmentation cache is required; legacy responses.h5 is not a cache')
        fit_settings = json.loads(source.attrs['effective_settings'])
        fit_config = json.loads(source.attrs['config'])
        # Direct callers must not silently reinterpret a cache using different fitted data.
        requested = fitting_config(config)
        for key in fit_config:
            if key not in config:
                continue
            expected = fit_settings.get(key, fit_config[key])
            value = requested[key]
            if key == 'window_s' and value is None:
                continue
            if json.dumps(value) != json.dumps(expected):
                raise ValueError(f'{key} differs from the segmentation cache; regenerate fits.h5')
        settings = Settings(**fit_settings, **timing_config(config))
        for path in (output_h5, output_csv):
            Path(path).parent.mkdir(parents=True, exist_ok=True)
        with h5py.File(output_h5, 'w') as out:
            out.attrs.update(dict(source.attrs))
            del out.attrs['artifact_type']
            del out.attrs['complete']
            out.attrs['schema_version'] = 1
            out.attrs['effective_settings'] = json.dumps(asdict(settings))
            out.attrs['config'] = json.dumps({**fit_config, **timing_config(config)})
            out.attrs['segmentation_implementation_sha256'] = source.attrs['implementation_sha256']
            out.attrs['implementation_sha256'] = json.dumps({name: hashlib.sha256(
                Path(__file__).with_name(name).read_bytes()).hexdigest()
                for name in ('session.py', 'bayesian_blocks.py', 'settings.py', 'metadata.py')})
            out.attrs['segmentation_cache'] = json.dumps({**_input_record(fits_path), 'sha256': cache_hash})
            out.attrs['selection_backend'] = 'cached_edges_cpu'
            for name in ('sound_events', 'event_exclusion_flags'):
                source.copy(name, out)
            groups = out.create_group('conditions')
            for condition in json.loads(source.attrs['condition_names']):
                cached_group = source['conditions'][condition]
                group = groups.create_group(condition)
                group.attrs.update(dict(cached_group.attrs))
                for name in ('candidate_event_indices', 'event_indices', 'event_times_s'):
                    cached_group.copy(name, group)
                units = group.create_group('units')
                for name in json.loads(source.attrs['unit_names']):
                    cached_unit = cached_group['units'][name]
                    attrs = dict(cached_unit.attrs)
                    raw = []
                    for index in sorted(cached_unit['fits'], key=int):
                        fit = cached_unit['fits'][index]
                        raw.append({**dict(fit.attrs), 'raw_edges_s': fit['raw_edges_s'][:]})
                    if attrs['status'] == 'fitted':
                        if [fit['p0'] for fit in raw] != list(settings.p0_values):
                            raise ValueError(f'{condition}/{name}: incomplete penalty cache; regenerate fits.h5')
                    elif raw:
                        raise ValueError(f'{condition}/{name}: unavailable case has unexpected fits')
                    result = select_timing(raw, attrs['n_events'], attrs['n_spikes'],
                                           float(group.attrs['duration_s']), settings,
                                           status=attrs['status'], diagnostic=attrs['diagnostic'])
                    fits = result.pop('fits')
                    unit = units.create_group(name)
                    unit.attrs.update(result)
                    a = cached_unit['aligned_spikes_s'][:]
                    if len(a) != attrs['n_spikes']:
                        raise ValueError(f'{condition}/{name}: cached spike count mismatch')
                    edges, counts, rate = diagnostic_psth(a, attrs['n_events'], settings)
                    for key, data in [('psth_edges_s', edges), ('psth_counts', counts), ('psth_rate_hz', rate)]:
                        unit.create_dataset(key, data=data)
                    fit_group = unit.create_group('fits')
                    for index, fit in enumerate(fits):
                        fg = fit_group.create_group(str(index))
                        fg.attrs['p0'] = fit['p0']
                        for key in ('raw_edges_s', 'edges_s'):
                            fg.create_dataset(key, data=fit[key])
                    summary.append(dict(condition=condition, condition_label=group.attrs['label'], unit=name,
                                        sound_ids=json.dumps(group.attrs['sound_ids'].tolist()),
                                        duration_s=group.attrs['duration_s'], duration_status=group.attrs['duration_status'],
                                        n_events_total=group.attrs['n_events_total'],
                                        n_events_excluded=group.attrs['n_events_excluded'], **result))
                print(f'{condition}: selected timing for {len(units)} units from cached fits', flush=True)
            out.attrs['selection_elapsed_seconds'] = time.perf_counter() - started
    columns = ['condition', 'condition_label', 'unit', 'sound_ids', 'duration_s', 'duration_status',
               'n_events_total', 'n_events_excluded'] + [key for key in _empty_result('', 0, 0) if key != 'fits']
    with open(output_csv, 'w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(summary)
    return summary


def run_session(meta_path, units_path, output_h5, output_csv, config):
    """Convenience API for both stages; Snakemake schedules the stages separately."""
    from .fits_session import fit_session
    cache = Path(output_h5).with_suffix('.fits.h5')
    fit_session(meta_path, units_path, cache, config)
    return select_session(cache, output_h5, output_csv, config)
