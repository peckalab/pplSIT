"""Coordinate-aware synchronization of keypoint-MoSeq headings."""
import numpy as np


def interpolate_heading(times, timestamps, heading, position):
    """Resample image-frame radians onto the position data's arena axes.

    Unwrap before interpolation to follow the shortest arc across +/- pi.
    As in the other synchronized columns, samples outside the video interval
    retain the nearest endpoint. Consecutive frames must differ by less than
    pi for the direction of rotation to be identifiable.

    Position conversion is (arena_origin - pixel) * scale * flip_sign.
    Its derivative is therefore +scale when flip=True, otherwise -scale.
    Apply those signs to the heading vector, then return atan2(y, x) radians.
    """
    times = np.asarray(times, dtype=float)
    timestamps = np.asarray(timestamps, dtype=float)
    heading = np.asarray(heading, dtype=float)
    if (times.ndim != 1 or timestamps.ndim != 1 or heading.ndim != 1
            or not len(timestamps) or len(timestamps) != len(heading)):
        raise ValueError('Heading requires matching, nonempty 1-D video samples.')
    if not all(np.all(np.isfinite(a)) for a in (times, timestamps, heading)):
        raise ValueError('Heading synchronization requires finite times and angles.')
    if np.any(np.diff(timestamps) <= 0):
        raise ValueError('Video timestamps must be strictly increasing.')

    interpolated = np.interp(times, timestamps, np.unwrap(heading))
    x_sign = 1 if position['flip_x'] else -1
    y_sign = 1 if position['flip_y'] else -1
    return np.arctan2(y_sign * np.sin(interpolated),
                      x_sign * np.cos(interpolated))
