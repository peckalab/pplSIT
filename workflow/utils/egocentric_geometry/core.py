"""Geometry only: no neural response selection or tuning statistics."""
import numpy as np


DEFAULTS = dict(
    references=['head', 'body'],
    origin_landmarks=['left_ear', 'right_ear'],
    axes=dict(head=dict(front=['nose'], back=['left_ear', 'right_ear']),
              body=dict(front=['neck'], back=['tail_base'])),
    angular_bins=12,
    ignored_channels=[8],
    speaker_positions_m=None,
)

# Flags can overlap. Invalid results are NaN, never a fictitious zero bearing.
FLAGS = {
    1: 'outside_tracking',
    8: 'nonfinite_pose', 16: 'zero_length_axis',
    32: 'unavailable_moseq_heading', 64: 'unresolved_speaker',
    128: 'at_speaker', 256: 'outside_metadata_timeline',
}


def settings(config):
    result = dict(DEFAULTS)
    unknown = set(config)-set(result)
    if unknown:
        raise ValueError(f'Unknown geometry settings: {sorted(unknown)}')
    result.update(config)
    refs = result['references']
    if (not isinstance(refs, list) or not refs or len(set(refs)) != len(refs)
            or not set(refs) <= {'head', 'body', 'moseq_heading'}):
        raise ValueError('references must select head, body and/or moseq_heading without duplicates')
    for names in [result['origin_landmarks']]+[
            result['axes'][ref][end] for ref in refs if ref != 'moseq_heading'
            for end in ('front', 'back')]:
        if (not isinstance(names, list) or not names
                or any(not isinstance(n, str) or not n for n in names)
                or len(set(names)) != len(names)):
            raise ValueError('Landmark lists must contain distinct nonempty names')
    n = result['angular_bins']
    if isinstance(n, bool) or not isinstance(n, int) or n < 2:
        raise ValueError('angular_bins must be an integer >=2')
    if (not isinstance(result['ignored_channels'], list)
            or any(isinstance(x, bool) or not isinstance(x, int)
                   for x in result['ignored_channels'])):
        raise ValueError('ignored_channels must contain integers')
    return result


def wrap(angle):
    return (angle+np.pi) % (2*np.pi)-np.pi


def sample_indices(query, timestamps):
    """Nearest existing synchronized sample; ties choose the earlier sample.

    This does not inspect native video confidence or gaps, or resample poses.
    Queries outside the synchronized recording are marked unavailable.
    """
    t = np.asarray(timestamps, dtype=float)
    q = np.asarray(query, dtype=float)
    if (t.ndim != 1 or not len(t) or not np.all(np.isfinite(t))
            or np.any(np.diff(t) <= 0) or q.ndim != 1 or not np.all(np.isfinite(q))):
        raise ValueError('Tracking times must be finite, 1-D and strictly increasing')
    hi = np.searchsorted(t, q).clip(0, len(t)-1)
    lo = np.maximum(hi-1, 0)
    indices = np.where(np.abs(q-t[lo]) <= np.abs(t[hi]-q), lo, hi)
    flags = np.zeros(len(q), dtype=np.uint16)
    flags[(q < t[0]) | (q > t[-1])] |= 1
    return indices, flags


def sample_pose(query, timestamps, coordinates):
    """Read synchronized positions without adding a tracking-quality policy."""
    indices, support = sample_indices(query, timestamps)
    points, flags = {}, {}
    for name, xy in coordinates.items():
        xy = np.asarray(xy, dtype=float)
        if xy.shape != (len(timestamps), 2):
            raise ValueError(f'Mismatched tracking arrays for {name}')
        f = support.copy()
        f[~np.isfinite(xy[indices]).all(axis=1)] |= 8
        points[name] = xy[indices]
        flags[name] = f
    return points, flags


def pose_reference(points, point_flags, reference, config, moseq_heading=None):
    def mean(names):
        return np.mean([points[n] for n in names], axis=0)

    names = list(config['origin_landmarks'])
    origin = mean(names)
    if reference == 'moseq_heading':
        heading, extra_flags = moseq_heading
        heading = np.asarray(heading).copy()
    else:
        spec = config['axes'][reference]
        names += spec['front']+spec['back']
        axis = mean(spec['front'])-mean(spec['back'])
        heading = np.arctan2(axis[:, 1], axis[:, 0])
        extra_flags = np.zeros(len(origin), dtype=np.uint16)
        extra_flags[np.linalg.norm(axis, axis=1) == 0] |= 16
    flags = extra_flags.copy()
    for name in set(names):
        flags |= point_flags[name]
    heading[flags != 0] = np.nan
    origin[flags != 0] = np.nan
    return origin, wrap(heading), flags


def resolve_speaker(spec, positions, ignored_channels):
    channels = spec.get('channels')
    if not isinstance(channels, (list, tuple)):
        return None, 'missing_channels', None
    if any(isinstance(ch, bool) or not isinstance(ch, int) for ch in channels):
        return None, 'invalid_channels', None
    physical = sorted(set(channels)-set(ignored_channels))
    if not physical:
        return None, 'no_physical_source', None
    if len(physical) != 1:
        return None, 'multiple_sources', None
    channel = physical[0]
    if channel not in positions:
        return None, 'unknown_speaker', channel
    return np.asarray(positions[channel], dtype=float), 'single_source', channel


def sound_bearing(origin, heading, pose_flags, speaker_xy):
    delta = speaker_xy-origin
    flags = pose_flags.copy()
    flags[~np.isfinite(speaker_xy).all(axis=1)] |= 64
    distance = np.linalg.norm(delta, axis=1)
    flags[distance == 0] |= 128
    angle = wrap(np.arctan2(delta[:, 1], delta[:, 0])-heading)
    angle[flags != 0] = np.nan
    distance[flags != 0] = np.nan
    return angle, distance, flags
