"""Geometry invariants using existing synchronized tracking; no tracking QC."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'workflow'))
from utils.egocentric_geometry.core import (
    sample_indices, sample_pose, pose_reference, resolve_speaker, settings, sound_bearing)
from utils.egocentric_geometry.session import run_session, load_moseq, sample_moseq
from utils.egocentric_geometry.plots import plot_geometry
from utils.setup import get_channel_positions


class GeometryTests(unittest.TestCase):
    def test_nearest_sample_ties_and_no_extrapolation(self):
        indices, flags = sample_indices([-.1, 0, .05, .1, .5, 1, 1.1], [0, .1, 1])
        np.testing.assert_array_equal(flags, [1, 0, 0, 0, 0, 0, 1])
        np.testing.assert_array_equal(indices, [0, 0, 0, 1, 1, 2, 2])

    def test_existing_coordinates_are_sampled_without_smoothing(self):
        points, flags = sample_pose([0, .001, .01], [0, .01],
            {'nose': np.array([[0., 0.], [1., 0.]])})
        np.testing.assert_array_equal(flags['nose'], [0, 0, 0])
        np.testing.assert_array_equal(points['nose'], [[0, 0], [0, 0], [1, 0]])

    def test_undefined_axis_is_flagged_without_minimum_length_filter(self):
        config = settings({})
        points = {'left_ear': np.zeros((2, 2)), 'right_ear': np.zeros((2, 2)),
                  'nose': np.array([[0., 0.], [1e-8, 0.]])}
        flags = {n: np.zeros(2, dtype='uint16') for n in points}
        _, heading, result = pose_reference(points, flags, 'head', config)
        self.assertEqual(result[0], 16)
        self.assertTrue(np.isnan(heading[0]))
        self.assertEqual(result[1], 0)
        self.assertEqual(heading[1], 0)

    def test_cardinal_bearings_and_unresolved_source(self):
        origin = np.zeros((5, 2))
        speakers = np.array([[1, 0], [0, 1], [-1, 0], [0, -1], [np.nan, np.nan]])
        angle, distance, flags = sound_bearing(origin, np.zeros(5), np.zeros(5, dtype='uint16'), speakers)
        np.testing.assert_allclose(angle[:4], [0, np.pi/2, -np.pi, -np.pi/2])
        np.testing.assert_array_equal(distance[:4], 1)
        self.assertTrue(np.isnan(angle[-1]))
        self.assertEqual(flags[-1], 64)
        _, _, flags = sound_bearing(origin[:1], np.zeros(1), np.zeros(1, dtype='uint16'), origin[:1])
        self.assertEqual(flags[0], 128)

    def test_speaker_resolution_never_picks_arbitrary_channel(self):
        positions = {3: [0, .5], 4: [.5, 0]}
        xy, status, channel = resolve_speaker({'channels': [8, 3]}, positions, [8])
        self.assertEqual(status, 'single_source')
        self.assertEqual(channel, 3)
        np.testing.assert_array_equal(xy, [0, .5])
        for spec, expected in [({}, 'missing_channels'), ({'channels': [8]}, 'no_physical_source'),
                               ({'channels': [3, 4, 8]}, 'multiple_sources'),
                               ({'channels': [9, 8]}, 'unknown_speaker')]:
            self.assertEqual(resolve_speaker(spec, positions, [8])[1], expected)

    def test_dates_and_copy_of_setup_positions(self):
        first = get_channel_positions('000001_SIT_2023-03-01_12-00-00')
        second = get_channel_positions('000001_SIT_2024-05-01_12-00-00')
        third = get_channel_positions('000001_SIT_2026-03-01_12-00-00')
        self.assertIn(1, first)
        self.assertLess(second[5][1], 0)
        self.assertGreater(third[5][1], 0)
        third[3] = [9, 9]
        self.assertEqual(get_channel_positions('000001_SIT_2026-03-01_12-00-00')[3], (0., .5))

    def test_modeled_heading_uses_existing_samples_without_interpolation(self):
        heading, flags = sample_moseq([.005, .5, 1.1], (np.array([0, .01, 1]),
                                      np.deg2rad([179, -179, -150])))
        np.testing.assert_allclose(heading[:2], np.deg2rad([179, -179]))
        np.testing.assert_array_equal(flags, [0, 0, 32])
        self.assertTrue(np.isnan(heading[2]))

    def test_invalid_settings(self):
        for config in [{'references': []}, {'references': ['head', 'head']},
                       {'references': ['unknown']}, {'likelihood_threshold': .9},
                       {'angular_bins': 2.5}, {'max_video_gap_s': .1},
                       {'ignored_channels': ['8']}, {'typo': True}]:
            with self.subTest(config=config), self.assertRaises(ValueError):
                settings(config)

    def test_session_io_active_and_passive(self):
        for passive in (False, True):
            with self.subTest(passive=passive), tempfile.TemporaryDirectory() as folder:
                root = Path(folder)
                position = dict(arena_x=100, arena_y=100, floor_r_in_meters=.5,
                                floor_radius=50, flip_x=True, flip_y=True)
                sounds = ({'F1000': {'channels': [3, 8]}, 'L2': {'channels': [3, 4, 8]}}
                          if passive else {'background': {'channels': [3, 8]},
                                           'target': {'channels': [3, 4, 8]}})
                parameters = dict(position=position, sound=dict(sounds=sounds))
                timeline = np.arange(0, .61, .01)
                coordinates = dict(left_ear=[0, .01], right_ear=[0, -.01],
                                   nose=[.02, 0], neck=[0, .08], tail_base=[0, -.02])
                headers = ['time']+[f'{part}_{c}' for part in coordinates for c in ['x', 'y', 'likelihood']]
                # Even zero confidence is retained: this stage adds no quality filter.
                data = np.column_stack([timeline, np.tile(np.concatenate(
                    [xy+[0.] for xy in coordinates.values()]), (len(timeline), 1))])
                data[10, headers.index('nose_x')] = np.nan
                data[30, headers.index('neck_y')] = -.02  # zero-length body axis
                with h5py.File(root/'dlc.h5', 'w') as f:
                    ds = f.create_dataset('df_DLC', data=data)
                    ds.attrs['headers'] = headers
                with h5py.File(root/'meta.h5', 'w') as f:
                    p = f.create_group('processed')
                    p.attrs['parameters'] = json.dumps(parameters)
                    p.create_dataset('timeline', data=timeline[:, None])
                    p.create_dataset('sound_events', data=[
                        [0, 1, 0], [.05, 1, 5], [.1, 1, 10], [.3, 1, 30],
                        [.5, 2, 50], [.7, 1, 70], [.55, 0, 55]])
                inputs = {key: str(root/name) for key, name in dict(meta='meta.h5', dlc='dlc.h5').items()}
                inputs['moseq'] = None
                rows = run_session(inputs, root/'geometry.h5', root/'summary.csv',
                    {}, '000001_SIT_2026-03-01_12-00-00')
                self.assertEqual(len(rows), 4)
                with h5py.File(root/'geometry.h5', 'r') as f:
                    self.assertTrue(f.attrs['complete'])
                    np.testing.assert_array_equal(f['event_indices'][:], np.arange(6))
                    self.assertAlmostEqual(f['events/head/bearing_rad'][0], np.pi/2)
                    self.assertAlmostEqual(f['events/body/bearing_rad'][0], 0)
                    np.testing.assert_array_equal(f['events/head/origin_xy_m'][0],
                                                  f['events/body/origin_xy_m'][0])
                    self.assertEqual(f['events/head/pose_flags'][1], 0)
                    self.assertTrue(f['events/head/pose_flags'][2] & 8)
                    self.assertEqual(f['events/body/pose_flags'][1], 0)
                    self.assertTrue(f['events/body/pose_flags'][3] & 16)
                    self.assertTrue(f['events/body/geometry_flags'][4] & 64)
                    self.assertTrue(f['events/body/geometry_flags'][5] & 256)
                    self.assertEqual(f['events/head/coverage/1/counts'][:].sum(), 3)
                    np.testing.assert_allclose(f['event_tracking_time_offsets_s'][:5], 0, atol=1e-12)
                plot_geometry(root/'geometry.h5', root/'coverage.pdf', root/'coverage.png')
                self.assertTrue((root/'coverage.pdf').read_bytes().startswith(b'%PDF'))

    def test_stale_moseq_is_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/'MoSeq.h5'
            with h5py.File(path, 'w') as f:
                f.create_dataset('moseq', data=np.zeros((2, 2)))
            with self.assertRaisesRegex(ValueError, 'regenerate'):
                load_moseq(path, {'flip_x': True, 'flip_y': False})


if __name__ == '__main__':
    unittest.main()
