"""MoSeq heading regressions; no model inference, GPU or Snakemake required."""
import itertools
import json
from pathlib import Path
import runpy
import sys
import tempfile
from types import SimpleNamespace
import unittest

import h5py
import numpy as np
import pandas as pd

WORKFLOW = Path(__file__).resolve().parents[1] / 'workflow'
sys.path.insert(0, str(WORKFLOW))
from utils.moseq import interpolate_heading


class HeadingTests(unittest.TestCase):
    def assert_angles(self, actual, expected):
        np.testing.assert_allclose(np.cos(actual), np.cos(expected), atol=1e-12)
        np.testing.assert_allclose(np.sin(actual), np.sin(expected), atol=1e-12)

    def test_wrap_crossings_both_directions(self):
        for degrees in ([170, -170], [-170, 170]):
            with self.subTest(degrees=degrees):
                actual = interpolate_heading([0, .5, 1], [0, 1],
                    np.deg2rad(degrees), dict(flip_x=True, flip_y=True))
                self.assert_angles(actual, np.deg2rad([degrees[0], 180, degrees[1]]))
                self.assertTrue(np.all(np.abs(actual) <= np.pi))

    def test_all_axis_flips_match_transformed_landmark_vectors(self):
        # Independent oracle: transform landmark positions exactly as the
        # synchronization script does, then measure their direction.
        front_x = np.array([12., 8., 8., 12.])
        front_y = np.array([21., 21., 19., 19.])
        back_x, back_y = 10., 20.
        heading = np.arctan2(front_y-back_y, front_x-back_x)
        t = np.arange(4.)
        for flip_x, flip_y in itertools.product([False, True], repeat=2):
            with self.subTest(flip_x=flip_x, flip_y=flip_y):
                x = (100-front_x)*.01*(-1 if flip_x else 1)
                y = (200-front_y)*.01*(-1 if flip_y else 1)
                bx = (100-back_x)*.01*(-1 if flip_x else 1)
                by = (200-back_y)*.01*(-1 if flip_y else 1)
                actual = interpolate_heading(t, t, heading,
                                             dict(flip_x=flip_x, flip_y=flip_y))
                self.assert_angles(actual, np.arctan2(y-by, x-bx))

    def test_pilot_flip_and_wrap_together(self):
        actual = interpolate_heading([0, .5, 1], [0, 1], np.deg2rad([170, -170]),
                                     dict(flip_x=True, flip_y=False))
        self.assert_angles(actual, np.deg2rad([-170, 180, 170]))

    def test_multiple_turns_and_irregular_timestamps(self):
        timestamps = [0, 1, 3, 4, 6]
        wrapped = np.deg2rad([0, 120, -120, 0, 120])
        actual = interpolate_heading([.5, 2, 3.5, 5], timestamps, wrapped,
                                     dict(flip_x=True, flip_y=True))
        self.assert_angles(actual, np.deg2rad([60, 180, 300, 420]))

    def test_endpoint_hold_and_single_frame(self):
        actual = interpolate_heading([-1, 2], [0, 1], [.2, .6],
                                     dict(flip_x=True, flip_y=False))
        self.assert_angles(actual, [-.2, -.6])
        actual = interpolate_heading([-1, 0, 2], [0], [.2],
                                     dict(flip_x=True, flip_y=False))
        self.assert_angles(actual, [-.2]*3)

    def test_invalid_samples_rejected(self):
        cases = [
            ([0], [], []), ([0], [0, 1], [0]),
            ([0], [0, 0], [0, 1]), ([0], [1, 0], [0, 1]),
            ([0], [0, np.nan], [0, 1]), ([0], [0, 1], [0, np.inf]),
            ([np.nan], [0, 1], [0, 1]), ([0], [0, 1], [[0, 1]]),
        ]
        for times, timestamps, angles in cases:
            with self.subTest(case=(times, timestamps, angles)), self.assertRaises(ValueError):
                interpolate_heading(times, timestamps, angles,
                                    dict(flip_x=True, flip_y=False))

    def test_actual_script_preserves_other_columns_and_records_convention(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            config = {'position': dict(arena_x=100, arena_y=200,
                floor_r_in_meters=.5, floor_radius=100, flip_x=True, flip_y=False),
                'video': {'fps': 50}}
            (root/'session.json').write_text(json.dumps(config))
            np.savetxt(root/'events.csv', [[100., 0], [100.05, 1]],
                       delimiter=',', header='time,event', comments='')
            np.savetxt(root/'raw.csv', [100., 100.02, 100.04])
            source = pd.DataFrame({
                'syllables reindexed': [1, 4, 2],
                'syllables non-reindexed': [7, 9, 8],
                'centroid x': [102., 104., 106.],
                'centroid y': [201., 202., 204.],
                'heading': np.deg2rad([170, -170, -150]),
                'estimated bodypart0 x': [105., 108., 107.],
                'estimated bodypart0 y': [204., 206., 205.],
                'latent_state 0': [1., .5, -1.],
            })
            source.to_csv(root/'moseq.csv', index=False)
            fake = SimpleNamespace(
                scriptdir=str(WORKFLOW/'scripts'),
                input=SimpleNamespace(session_cfg=root/'session.json',
                    events=root/'events.csv', timestamps_raw=root/'raw.csv',
                    moseq_csv_path=root/'moseq.csv'),
                output=SimpleNamespace(moseq_h5_path_100Hz=root/'MoSeq.h5'))
            runpy.run_path(str(WORKFLOW/'scripts/synch_moseq_100Hz.py'),
                           init_globals={'snakemake': fake})
            with h5py.File(root/'MoSeq.h5', 'r') as f:
                ds = f['moseq']
                headers = list(ds.attrs['headers'])
                actual = ds[:]
                self.assertEqual(headers, ['time']+list(source.columns))
                self.assertEqual(ds.attrs['heading_processing'], 'arena_circular_v1')
                self.assertEqual(ds.attrs['heading_coordinate_system'], 'arena')
                self.assertEqual(ds.attrs['heading_units'], 'radians')
                self.assertTrue(ds.attrs['heading_flip_x'])
                self.assertFalse(ds.attrs['heading_flip_y'])
            times = np.arange(0, 100.05-100., .01)
            stamps = np.loadtxt(root/'raw.csv')-100.
            np.testing.assert_array_equal(actual[:, 0], times)
            nearest = [np.abs(stamps-t).argmin() for t in times]
            for col in source:
                if col == 'heading':
                    self.assert_angles(actual[:, headers.index(col)],
                                       np.deg2rad([-170, -180, 170, 160, 150]))
                    continue
                if 'syllables' in col:
                    expected = source[col].to_numpy()[nearest]
                else:
                    expected = np.interp(times, stamps, source[col])
                if ' x' in col:
                    expected = (100-expected)*.005*(-1)
                elif ' y' in col:
                    expected = (200-expected)*.005
                np.testing.assert_array_equal(actual[:, headers.index(col)], expected)


if __name__ == '__main__':
    unittest.main()
