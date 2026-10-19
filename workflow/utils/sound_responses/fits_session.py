"""Session I/O for the persistent, timing-independent segmentation cache."""
from dataclasses import asdict
import hashlib
import importlib.metadata
import json
from pathlib import Path
import re

import h5py
import numpy as np

from utils.session import build_stimulus_catalog
from .settings import Settings, FIT_SETTING_KEYS, fitting_config
from .segmentation import align_spikes, fit_spikes, select_events
from .metadata import METHOD, METHOD_VERSION, REFERENCE, CACHE_ARTIFACT, CACHE_SCHEMA_VERSION


def make_conditions(events, catalog, configured=None, include_pooled=True):
    """Positive stimulus IDs represent conditions; channel numbers are metadata."""
    observed = sorted({int(x) for x in events[events[:, 1] > 0, 1]})
    if configured is None:
        mapping = {f'sound_{sid}': [sid] for sid in observed}
    else:
        if not isinstance(configured, dict) or not configured:
            raise ValueError('conditions must be null or a nonempty mapping from names to sound IDs')
        mapping = dict(configured)
    # Validate custom IDs before constructing a pooled condition.
    for name, ids in mapping.items():
        if (not isinstance(ids, (list, tuple)) or not ids
                or any(isinstance(x, bool) or not isinstance(x, int) or x <= 0 for x in ids)):
            raise ValueError(f'{name}: sound IDs must be a nonempty list of positive integers')
    if include_pooled and observed:
        if 'pooled' in mapping:
            raise ValueError('pooled is reserved when include_pooled is true')
        # Pool the explicitly selected stimuli if custom conditions are given.
        pooled = sorted(set(sid for ids in mapping.values() for sid in ids))
        mapping['pooled'] = pooled
    conditions = []
    for name, ids in mapping.items():
        if not isinstance(name, str) or not re.fullmatch(r'[A-Za-z0-9_-]+', name):
            raise ValueError('Condition names must contain only letters, digits, underscores or hyphens')
        if len(set(ids)) != len(ids):
            raise ValueError(f'{name}: duplicate sound IDs')
        missing = set(ids) - set(observed)
        if missing:
            raise ValueError(f'{name}: sound IDs absent from session: {sorted(missing)}')
        specs = [catalog.get(str(sid), {}) for sid in ids]
        durations = [s.get('duration') for s in specs]
        duration = np.nan
        duration_status = 'missing_duration'
        if all(isinstance(x, (int, float)) and np.isfinite(x) and x > 0 for x in durations):
            duration_status = 'mixed_durations'
            if np.allclose(durations, durations[0], rtol=0, atol=1e-12):
                duration = float(durations[0])
                duration_status = 'known'
        label = specs[0].get('name', name) if configured is None and len(ids) == 1 and name != 'pooled' else name
        conditions.append(dict(name=name, label=label, sound_ids=ids, duration_s=duration,
                               duration_status=duration_status, stimuli=specs))
    return conditions


def _input_record(path):
    path = Path(path).resolve()
    stat = path.stat()
    return dict(path=str(path), size_bytes=stat.st_size, mtime_ns=stat.st_mtime_ns)


def fit_session(meta_path, units_path, output_h5, config):
    """Persist all raw penalty levels and aligned spikes for independent selection."""
    config = fitting_config(config)
    with h5py.File(meta_path, 'r') as meta:
        parameters = json.loads(meta['processed'].attrs['parameters'])
        catalog_json = meta.attrs.get('stimulus_catalog', meta['processed'].attrs.get('stimulus_catalog'))
        catalog = json.loads(catalog_json) if catalog_json is not None else build_stimulus_catalog(parameters)
        events = meta['processed/sound_events'][:]
        timeline = meta['processed/timeline']
        if len(timeline) < 2:
            raise ValueError('meta timeline does not establish recording bounds')
        bounds = (float(timeline[0, 0]), float(timeline[-1, 0]))
        sync_mode = str(meta.attrs.get('sounds_sync_mode', 'unrecorded'))
    settings = Settings.from_config(config, parameters)
    backend_name = config.get('backend', 'cpu')
    if backend_name not in ('cpu', 'gpu'):
        raise ValueError('sound_responses_bayesian_blocks.backend must be cpu or gpu')
    gpu = None
    if backend_name == 'gpu':
        from .gpu import GPUBackend
        gpu = GPUBackend(settings, memory_mb=config.get('gpu_memory_mb', 4096))
    eligible, exclusion = select_events(events, settings.window_s, bounds, config.get('time_ranges_s'))
    conditions = make_conditions(events, catalog, config.get('conditions'), config.get('include_pooled', True))
    versions = {name: importlib.metadata.version(name) for name in ('hepstats', 'numpy', 'pandas', 'h5py')}
    implementation_files = ['settings.py', 'segmentation.py', 'fits_session.py', 'metadata.py']
    if gpu is not None:
        implementation_files += ['gpu.py', 'gpu_kernel.cu']
    implementation = {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                      for name in implementation_files}
    Path(output_h5).parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(units_path, 'r') as units, h5py.File(output_h5, 'w') as out:
        unit_names = config.get('unit_names')
        if unit_names is None:
            unit_names = sorted(name for name in units if isinstance(units[name], h5py.Group) and 'spike_times' in units[name])
        if not isinstance(unit_names, list) or len(set(unit_names)) != len(unit_names):
            raise ValueError('unit_names must be null or a list without duplicates')
        for name in unit_names:
            if not isinstance(name, str) or name not in units or 'spike_times' not in units[name]:
                raise ValueError(f'Unknown unit: {name}')
        # Across-condition batching keeps the GPU occupied even for smaller populations.
        # These aligned spikes are also reused for the diagnostic PSTHs below.
        aligned, gpu_results = {}, {}
        if gpu is not None:
            cases, keys = [], []
            for condition in conditions:
                selected = np.flatnonzero(np.isin(events[:, 1], condition['sound_ids']) & eligible)
                for name in unit_names:
                    key = (condition['name'], name)
                    a = align_spikes(units[name]['spike_times'][:], events[selected, 0], settings.window_s)
                    aligned[key] = a
                    keys.append(key)
                    cases.append((a, len(selected)))
            gpu_results = dict(zip(keys, gpu.fit(cases)))
        out.attrs.update(schema_version=CACHE_SCHEMA_VERSION, artifact_type=CACHE_ARTIFACT, complete=False, method=METHOD, method_version=METHOD_VERSION,
                         reference=json.dumps(REFERENCE), package_versions=json.dumps(versions),
                         implementation_sha256=json.dumps(implementation),
                         backend=backend_name,
                         backend_provenance=json.dumps(gpu.provenance if gpu is not None else
                                                       {'name': 'hepstats', 'version': versions['hepstats']}),
                         config=json.dumps(config), effective_settings=json.dumps({key: asdict(settings)[key] for key in FIT_SETTING_KEYS}),
                         input_files=json.dumps([_input_record(meta_path), _input_record(units_path)]),
                         stimulus_catalog=json.dumps(catalog), recording_bounds_s=bounds,
                         sounds_sync_mode=sync_mode, unit_names=json.dumps(unit_names),
                         condition_names=json.dumps([condition['name'] for condition in conditions]),
                         is_modulated_encoding='-1 unavailable; 0 no accepted change; 1 accepted change',
                         latency_reference='seconds from meta sound_events onset for BOTH onset and offset',
                         p0_interpretation='segmentation regularization; not a responsiveness-test p-value',
                         exposure_policy='complete fixed windows; exclude overlaps with next logged event',
                         interval_convention='[onset, onset + window_s)',
                         exclusion_flags='1 outside recording; 2 overlaps next event; 4 outside time_ranges_s')
        out.create_dataset('sound_events', data=events)
        out.create_dataset('event_exclusion_flags', data=exclusion)
        groups = out.create_group('conditions')
        for condition in conditions:
            group = groups.create_group(condition['name'])
            ids = condition['sound_ids']
            candidates = np.flatnonzero(np.isin(events[:, 1], ids))
            selected = candidates[eligible[candidates]]
            group.attrs.update(label=condition['label'], sound_ids=ids,
                               duration_s=condition['duration_s'], duration_status=condition['duration_status'],
                               stimuli=json.dumps(condition['stimuli']), n_events_total=len(candidates),
                               n_events_used=len(selected), n_events_excluded=len(candidates) - len(selected))
            group.create_dataset('candidate_event_indices', data=candidates)
            group.create_dataset('event_indices', data=selected)
            group.create_dataset('event_times_s', data=events[selected, 0])
            result_group = group.create_group('units')
            print(f"{condition['name']}: {len(selected)}/{len(candidates)} complete event windows", flush=True)
            for name in unit_names:
                key = (condition['name'], name)
                a = aligned.pop(key) if gpu is not None else align_spikes(
                    units[name]['spike_times'][:], events[selected, 0], settings.window_s)
                result = gpu_results.pop(key) if gpu is not None else fit_spikes(
                    a, len(selected), settings)
                fits = result.pop('fits')
                unit = result_group.create_group(name)
                unit.attrs.update(result)
                unit.create_dataset('aligned_spikes_s', data=a, compression='gzip', shuffle=True)
                fit_group = unit.create_group('fits')
                for i, fit in enumerate(fits):
                    fg = fit_group.create_group(str(i))
                    fg.attrs.update(p0=fit['p0'], status=fit['status'], diagnostic=fit['diagnostic'])
                    fg.create_dataset('raw_edges_s', data=fit['raw_edges_s'])
                print(f"  {name}: {result['status']}, raw fits={len(fits)}", flush=True)
        out.attrs['n_conditions'] = len(conditions)
        out.attrs['n_units'] = len(unit_names)
        out.attrs['complete'] = True
