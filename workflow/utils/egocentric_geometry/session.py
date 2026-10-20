"""Build reusable geometry from the pipeline's synchronized tracking outputs."""
import csv
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np

from utils.session import build_stimulus_catalog
from utils.setup import get_channel_positions
from .core import (FLAGS, sample_indices, sample_pose, pose_reference,
                   resolve_speaker, settings, sound_bearing, wrap)


def file_record(path):
    p = Path(path).resolve()
    with p.open('rb') as stream:
        digest = hashlib.file_digest(stream, 'sha256').hexdigest()
    return dict(path=str(p), sha256=digest, size_bytes=p.stat().st_size)


def load_dlc(path, config):
    needed = set(config['origin_landmarks'])
    for ref in config['references']:
        if ref != 'moseq_heading':
            needed.update(config['axes'][ref]['front']+config['axes'][ref]['back'])
    points = {}
    with h5py.File(path, 'r') as f:
        ds = f['df_DLC']
        headers = list(ds.attrs['headers'])
        timestamps = ds[:, headers.index('time')]
        for name in sorted(needed):
            for axis in ('x', 'y'):
                if name+'_'+axis not in headers:
                    raise ValueError(f'Missing synchronized DLC landmark/component: {name}/{axis}')
            points[name] = np.column_stack([ds[:, headers.index(name+'_'+axis)] for axis in ('x', 'y')])
    return timestamps, points


def load_moseq(path, position):
    with h5py.File(path, 'r') as f:
        ds = f['moseq']
        if ds.attrs.get('heading_processing') != 'arena_circular_v1':
            raise ValueError('MoSeq heading is not corrected: regenerate synch_moseq_100Hz first')
        for axis in ('x', 'y'):
            if ds.attrs.get('heading_flip_'+axis) != position['flip_'+axis]:
                raise ValueError('MoSeq heading calibration differs from session position calibration')
        headers = list(ds.attrs['headers'])
        return ds[:, headers.index('time')], ds[:, headers.index('heading')]


def sample_moseq(query, samples):
    times, heading = samples
    indices, support = sample_indices(query, times)
    flags = np.zeros(len(query), dtype=np.uint16)
    flags[(support != 0) | ~np.isfinite(heading[indices])] |= 32
    angle = wrap(heading[indices])
    angle[flags != 0] = np.nan
    return angle, flags


def run_session(inputs, output_h5, output_csv, config, session):
    config = settings(config)
    paths = [Path(p).resolve() for p in inputs.values() if p]
    if (Path(output_h5).resolve() in paths or Path(output_csv).resolve() in paths
            or Path(output_h5).resolve() == Path(output_csv).resolve()):
        raise ValueError('Geometry outputs must be distinct from inputs and each other')
    with h5py.File(inputs['meta'], 'r') as meta:
        all_events = meta['processed/sound_events'][:]
        timeline = meta['processed/timeline'][:, 0]
        parameters = json.loads(meta['processed'].attrs['parameters'])
        encoded = meta.attrs.get('stimulus_catalog', meta['processed'].attrs.get('stimulus_catalog'))
        catalog = json.loads(encoded) if encoded is not None else build_stimulus_catalog(parameters)
        sync_mode = str(meta.attrs.get('sounds_sync_mode', 'unrecorded'))
    if (all_events.ndim != 2 or all_events.shape[1] < 2
            or not np.all(np.isfinite(all_events[:, :2]))
            or np.any(all_events[:, 1] != np.floor(all_events[:, 1]))):
        raise ValueError('Sound events require finite onset times and integer stimulus IDs')
    if len(timeline) == 0 or np.any(np.diff(timeline) <= 0):
        raise ValueError('Metadata timeline must be nonempty and strictly increasing')
    event_indices = np.flatnonzero(all_events[:, 1] > 0)
    events = all_events[event_indices]
    timestamps, points = load_dlc(inputs['dlc'], config)
    sample_indices(timeline, timestamps)
    moseq = (load_moseq(inputs['moseq'], parameters['position'])
             if 'moseq_heading' in config['references'] else None)
    if moseq is not None and not np.array_equal(moseq[0], timestamps):
        raise ValueError('Synchronized DLC and MoSeq time grids differ; synchronize them consistently first')
    positions = get_channel_positions(session)
    if config['speaker_positions_m'] is not None:
        positions.update({int(k): v for k, v in config['speaker_positions_m'].items()})
    for channel, xy in positions.items():
        if np.asarray(xy).shape != (2,) or not np.all(np.isfinite(xy)):
            raise ValueError(f'Speaker {channel} needs two finite arena coordinates in meters')
    speakers = np.full((len(events), 2), np.nan)
    conditions = []
    for sound_id in sorted(set(events[:, 1].astype(int))):
        spec = catalog.get(str(sound_id), {})
        xy, status, channel = resolve_speaker(spec, positions, config['ignored_channels'])
        selected = events[:, 1] == sound_id
        if xy is not None:
            speakers[selected] = xy
        conditions.append(dict(sound_id=int(sound_id), name=spec.get('name', f'sound_{sound_id}'),
                               speaker_status=status, channel=channel, stimulus=spec))
    input_records = [file_record(p) for p in inputs.values() if p]
    source_files = [Path(__file__), Path(__file__).with_name('core.py'),
                    Path(__file__).parents[1]/'setup.py', Path(__file__).parents[1]/'session.py']
    implementation = {str(p.relative_to(Path(__file__).parents[2])): file_record(p)['sha256']
                      for p in source_files}
    edges = np.linspace(-np.pi, np.pi, config['angular_bins']+1)
    rows = []
    Path(output_h5).parent.mkdir(parents=True, exist_ok=True)
    Path(output_csv).parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(output_h5, 'w') as out:
        out.attrs.update(artifact_type='egocentric_sound_geometry', schema_version=1, complete=False,
                         session=session, config=json.dumps(config), references=json.dumps(config['references']),
                         input_files=json.dumps(input_records), implementation_sha256=json.dumps(implementation),
                         stimulus_catalog=json.dumps(catalog), conditions=json.dumps(conditions),
                         speaker_positions_m=json.dumps(positions), speaker_position_source='setup.py + configured overrides',
                         position_units='meters', angle_units='radians',
                         angle_convention='[-pi, pi); 0=in front; positive=counterclockwise in arena x/y',
                         flags=json.dumps(FLAGS), sounds_sync_mode=sync_mode,
                         time_reference='seconds since first behavioral event; sounds from meta sound_events',
                         pose_sampling='nearest existing synchronized DLC sample; ties choose earlier sample',
                         tracking_quality_policy='unchanged upstream tracking; no confidence filtering or native-video gap checks',
                         event_selection='all positive sound IDs; no neural response or Bayesian-block eligibility filtering')
        out.create_dataset('source_sound_events', data=all_events)
        out.create_dataset('event_indices', data=event_indices)
        out.create_dataset('event_times_s', data=events[:, 0])
        out.create_dataset('sound_ids', data=events[:, 1].astype(int))
        out.create_dataset('speaker_xy_m', data=speakers)
        out.create_dataset('timeline_s', data=timeline)
        event_samples, event_support = sample_indices(events[:, 0], timestamps)
        sample_times = timestamps[event_samples].copy()
        sample_times[event_support != 0] = np.nan
        out.create_dataset('event_tracking_times_s', data=sample_times)
        out.create_dataset('event_tracking_time_offsets_s', data=sample_times-events[:, 0])
        out.create_dataset('angular_bin_edges_rad', data=edges)
        timeline_results = {}
        for scope, query in [('timeline', timeline), ('events', events[:, 0])]:
            sampled, point_flags = sample_pose(query, timestamps, points)
            sample_idx, support = sample_indices(query, timestamps)
            sampled_times = timestamps[sample_idx].copy()
            sampled_times[support != 0] = np.nan
            out.require_group(scope).create_dataset('tracking_times_s', data=sampled_times)
            modeled = sample_moseq(query, moseq) if moseq is not None else None
            for ref in config['references']:
                origin, heading, flags = pose_reference(sampled, point_flags, ref, config, modeled)
                flags[(query < timeline[0]) | (query > timeline[-1])] |= 256
                heading[flags != 0] = np.nan
                origin[flags != 0] = np.nan
                group = out.require_group(scope).create_group(ref)
                group.create_dataset('origin_xy_m', data=origin, compression='gzip')
                group.create_dataset('heading_rad', data=heading, compression='gzip')
                group.create_dataset('pose_flags', data=flags, compression='gzip')
                if scope == 'timeline':
                    timeline_results[ref] = heading
                    continue
                angle, distance, geometry_flags = sound_bearing(origin, heading, flags, speakers)
                group.create_dataset('bearing_rad', data=angle)
                group.create_dataset('distance_m', data=distance)
                group.create_dataset('geometry_flags', data=geometry_flags)
                for condition in conditions:
                    subset = events[:, 1] == condition['sound_id']
                    valid = subset & (geometry_flags == 0)
                    histogram = np.histogram(angle[valid], bins=edges)[0]
                    coverage = group.require_group('coverage').create_group(str(condition['sound_id']))
                    coverage.create_dataset('counts', data=histogram)
                    coverage.attrs.update(n_events=int(subset.sum()), n_valid=int(valid.sum()),
                                          label=condition['name'], speaker_status=condition['speaker_status'])
                    row = dict(reference=ref, sound_id=condition['sound_id'], condition=condition['name'],
                               speaker_status=condition['speaker_status'], n_events=int(subset.sum()),
                               n_valid=int(valid.sum()), valid_fraction=float(valid.sum()/subset.sum()),
                               occupied_bins=int((histogram > 0).sum()), min_events_per_bin=int(histogram.min()))
                    row.update({name: int(((geometry_flags[subset] & bit) != 0).sum())
                                for bit, name in FLAGS.items()})
                    rows.append(row)
        comparison = {}
        refs = config['references']
        for i, a in enumerate(refs):
            for b in refs[i+1:]:
                valid = np.isfinite(timeline_results[a]) & np.isfinite(timeline_results[b])
                delta = np.abs(wrap(timeline_results[a][valid]-timeline_results[b][valid]))
                comparison[a+'__'+b] = dict(n_samples=int(valid.sum()),
                    median_abs_difference_deg=float(np.rad2deg(np.median(delta))) if len(delta) else None)
        out.attrs['reference_comparisons'] = json.dumps(comparison)
        out.attrs['complete'] = True
    fields = ['reference', 'sound_id', 'condition', 'speaker_status', 'n_events', 'n_valid',
              'valid_fraction', 'occupied_bins', 'min_events_per_bin']+list(FLAGS.values())
    with open(output_csv, 'w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    return rows
