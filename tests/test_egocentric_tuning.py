"""Reference fixtures, statistical invariants, and end-to-end event joins."""
import itertools
import json
from pathlib import Path
import sys
import tempfile
import unittest

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'workflow'))
from utils.egocentric_tuning.core import (
    CENTERS, PHASES, analyze, angle_bins, bootstrap_curves, compare, count_spikes,
    criteria, phase_windows, settings, stable_seed, vector_metrics)
from utils.egocentric_tuning.session import run_session
from utils.egocentric_tuning.plots import plot_tuning


class TuningTests(unittest.TestCase):
    def test_reference_observed_curves_and_metrics(self):
        fixture = json.loads(Path(__file__).with_name('fixtures').joinpath('amaro_centered_tuning.json').read_text())
        for case in fixture['cases']:
            result, _ = analyze(np.deg2rad(case['angles_deg']), case['counts'], .05, 20, 3)
            np.testing.assert_allclose(result['spikes_per_presentation'], case['curve'], atol=1e-14)
            np.testing.assert_allclose([result['vector_strength'], result['orientation_strength']],
                                       case['metrics'][:2], atol=1e-14)
            if result['vector_strength'] > 1e-12:
                self.assertAlmostEqual(np.rad2deg(result['preferred_direction_rad']), case['metrics'][2], places=12)
            if result['orientation_strength'] > 1e-12:
                self.assertAlmostEqual(np.rad2deg(result['preferred_orientation_rad']), case['metrics'][3], places=12)

    def test_centered_bins_and_occupancy_correction(self):
        np.testing.assert_array_equal(angle_bins(np.deg2rad([-22.5, 22.499, 22.5, 337.5, 360])), [0, 0, 1, 0, 0])
        angles = np.repeat(CENTERS, np.arange(10, 90, 10))
        r, _ = analyze(angles, np.full(len(angles), 3), .05, 50, 2)
        np.testing.assert_allclose(r['rate_hz'], 60)
        self.assertLess(r['vector_strength'], 1e-12)
        self.assertEqual(r['classification'], 'none')
        self.assertTrue(np.isnan(r['preferred_direction_rad']))

    def test_intervals_partition_observation_and_half_open_boundaries(self):
        windows = phase_windows(.01, .07, .05, .24)
        spikes = np.array([0., .009, .01, .069, .07, .119, .07+.05, .239, .24])
        expected = dict(onset=2, offset=2, late=4)
        for phase in PHASES:
            self.assertEqual(count_spikes(spikes, np.array([0.]), windows[phase][0])[0], expected[phase])
        self.assertAlmostEqual(sum(np.diff(w[0], axis=1).sum() for w in windows.values()), .24)
        fallback = phase_windows(.01, np.nan, .08, .24)
        np.testing.assert_allclose(fallback['onset'][0], [[.01, .09]])
        np.testing.assert_allclose(fallback['offset'][0], [[.09, .17]])
        self.assertEqual(fallback['offset'][1], 'duration_fallback')
        self.assertTrue(all(len(w[0]) == 0 for w in phase_windows(np.nan, .07, .05, .24).values()))
        self.assertTrue(all(len(w[0]) == 0 for w in phase_windows(.01, .07, np.nan, .24).values()))
        outside = phase_windows(.1, .22, .05, .24)
        self.assertEqual(len(outside['onset'][0]), 1)
        self.assertEqual(len(outside['offset'][0]), 0)
        self.assertEqual(len(outside['late'][0]), 0)

    def test_bootstrap_matches_exact_whole_event_resampling_distribution(self):
        # All 4^4 possible resamples; zero-spike events must stay in the pool.
        bins = np.array([0, 0, 0, 1])
        counts = np.array([0, 0, 3, 4])
        exact = []
        for draw in itertools.product(range(4), repeat=4):
            b, c = bins[list(draw)], counts[list(draw)]
            n = np.bincount(b, minlength=8)
            s = np.bincount(b, weights=c, minlength=8)
            exact.append(np.divide(s, n, out=np.zeros(8), where=n > 0))
        draws = bootstrap_curves(bins, counts, 100000, 11)
        np.testing.assert_allclose(draws.mean(axis=0), np.mean(exact, axis=0), atol=.015)
        np.testing.assert_allclose(draws.var(axis=0), np.var(exact, axis=0), atol=.025)
        np.testing.assert_array_equal(draws, bootstrap_curves(bins, counts, 100000, 11))
        self.assertNotEqual(stable_seed(0, 'a'), stable_seed(0, 'b'))

    def test_reference_criteria_and_strict_table_lookup(self):
        occupied = np.full(8, 10)
        r = criteria(75, occupied, [.2, .1], [.15, .1])
        self.assertEqual(r['classification'], 'directional')
        # Larger folded metric fails; smaller directional metric passing is not enough.
        r = criteria(100, occupied, [.3, .4], [.2, .1])
        self.assertTrue(r['directional_pass'])
        self.assertEqual(r['classification'], 'none')
        # Strict < row lookup at exactly 25 spikes takes the 20-spike row.
        r = criteria(25, occupied, [.9, .8], [.37, .4])
        self.assertAlmostEqual(r['directional_lower_threshold'], .3817631849352243)
        self.assertFalse(r['directional_pass'])
        self.assertEqual(criteria(20, occupied, [1, 1], [1, 1])['status'], 'insufficient_spikes')
        self.assertEqual(criteria(100, [9]*8, [1, 1], [1, 1])['classification'], 'unavailable')
        self.assertEqual(criteria(100, [101]+[10]*7, [1, 1], [1, 1])['status'], 'asymmetric_angular_coverage')

    def test_empty_and_zero_spike_events(self):
        r, _ = analyze([], [], .05, 20, 0)
        self.assertEqual(r['n_events'], 0)
        self.assertEqual(r['classification'], 'unavailable')
        self.assertTrue(np.all(np.isnan(r['rate_hz'])))
        r, _ = analyze(np.tile(CENTERS, 10), np.zeros(80), .05, 20, 0)
        self.assertEqual(r['n_events'], 80)
        self.assertEqual(r['status'], 'insufficient_spikes')
        np.testing.assert_array_equal(r['rate_hz'], 0)

    def test_comparison_sign_and_circular_shifts(self):
        angles = np.tile(CENTERS, 10)
        a, ba = analyze(angles, np.tile([8, 4, 0, 0, 0, 0, 0, 4], 10), .05, 100, 1)
        b, bb = analyze(angles, np.tile([4, 8, 4, 0, 0, 0, 0, 0], 10), .05, 100, 2)
        diff = compare(a, b, ba, bb)
        np.testing.assert_allclose(diff['rate_difference_hz'], b['rate_hz']-a['rate_hz'])
        self.assertAlmostEqual(diff['direction_shift_rad'], np.pi/4)
        self.assertAlmostEqual(diff['equal_angle_mean_difference_hz'], 0.)

    def test_config_validation(self):
        for config in ({'angular_bins': 12}, {'seed': -1}, {'references': 'head'},
                       {'n_bootstrap': True}, {'comparisons': [['a', 'a']]}):
            with self.assertRaises(ValueError):
                settings(config)

    def create_inputs(self, root):
        events = np.column_stack([np.arange(160)*.3, np.repeat([1, 2], 80), np.zeros(160)])
        units = root/'units.h5'
        with h5py.File(units, 'w') as f:
            # Same rate across angles; different frequency/task condition doubles counts.
            f.create_dataset('u/spike_times', data=np.sort(np.r_[events[:, 0]+.02, events[80:, 0]+.03]))
        stat = units.stat()
        with h5py.File(root/'responses.h5', 'w') as f:
            f.attrs.update(method='amaro_07062021_bayesian_blocks',
                           input_files=json.dumps([dict(path=str(units), size_bytes=stat.st_size, mtime_ns=stat.st_mtime_ns)]),
                           condition_names=json.dumps(['sound_1', 'sound_2', 'pooled']), unit_names='["u"]',
                           effective_settings='{"window_s": 0.24}', stimulus_catalog='{}')
            f.create_dataset('sound_events', data=events)
            for condition, indices in [('sound_1', np.arange(80)), ('sound_2', np.arange(80, 160)),
                                       ('pooled', np.arange(160))]:
                c = f.create_group('conditions/'+condition)
                c.attrs.update(duration_s=.05, duration_status='known', label=condition)
                c.create_dataset('event_indices', data=indices)
                c.create_dataset('event_times_s', data=events[indices, 0])
                u = c.create_group('units/u')
                u.attrs.update(onset_status='detected', onset_latency_s=.01, offset_status='detected',
                               offset_latency_s=.06, n_spikes=80 if condition == 'sound_1' else
                               (160 if condition == 'sound_2' else 240))
        with h5py.File(root/'geometry.h5', 'w') as f:
            f.attrs.update(artifact_type='egocentric_sound_geometry', complete=True,
                           references='["head"]', session='test')
            f.create_dataset('source_sound_events', data=events)
            f.create_dataset('event_indices', data=np.arange(160))
            f.create_dataset('event_times_s', data=events[:, 0])
            f.create_dataset('sound_ids', data=events[:, 1])
            f.create_dataset('events/head/bearing_rad', data=np.tile(CENTERS, 20))
            f.create_dataset('events/head/geometry_flags', data=np.zeros(160, dtype='uint16'))
        return {name: str(root/f'{name}.h5') for name in ('geometry', 'responses', 'units')}

    def test_end_to_end_and_plot(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            inputs = self.create_inputs(root)
            run_session(inputs, root/'tuning.h5', root/'summary.csv', root/'comparisons.csv',
                        dict(n_bootstrap=30), 'test')
            with h5py.File(root/'tuning.h5') as f:
                self.assertTrue(f.attrs['complete'])
                np.testing.assert_allclose(f['comparisons/0/u/head/onset/rate_difference_hz'][:], 20)
                self.assertEqual(f['conditions/sound_1/units/u/onset/head'].attrs['n_events'], 80)
            plot_tuning(root/'tuning.h5', root/'tuning.pdf', root/'preview.png')
            self.assertGreater((root/'tuning.pdf').stat().st_size, 1000)
            # Changed event order cannot silently join.
            with h5py.File(inputs['geometry'], 'a') as f:
                f['source_sound_events'][0, 0] = 1.
            with self.assertRaisesRegex(ValueError, 'sound events differ'):
                run_session(inputs, root/'bad.h5', root/'bad.csv', root/'bad2.csv', {}, 'test')

    def test_provenance_rejects_different_units(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            inputs = self.create_inputs(root)
            with h5py.File(inputs['units'], 'a') as f:
                f.attrs['changed'] = True
            with self.assertRaisesRegex(ValueError, 'Units export differs'):
                run_session(inputs, root/'bad.h5', root/'bad.csv', root/'bad2.csv', {}, 'test')

    def test_missing_timing_and_invalid_geometry_remain_explicit(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            inputs = self.create_inputs(root)
            with h5py.File(inputs['responses'], 'a') as f:
                f['conditions/sound_2/units/u'].attrs['onset_status'] = 'not_detected'
            with h5py.File(inputs['geometry'], 'a') as f:
                f['events/head/geometry_flags'][0] = 16
                f['events/head/bearing_rad'][0] = np.nan
            run_session(inputs, root/'tuning.h5', root/'summary.csv', root/'comparisons.csv',
                        dict(n_bootstrap=20, timing_condition=None), 'test')
            with h5py.File(root/'tuning.h5') as f:
                a = f['conditions/sound_1/units/u/onset/head']
                self.assertEqual(a.attrs['n_events'], 79)
                self.assertEqual(a.attrs['n_spikes'], 79)
                self.assertEqual(a.attrs['status'], 'insufficient_angular_coverage')
                b = f['conditions/sound_2/units/u/onset']
                self.assertEqual(b['head'].attrs['status'], 'timing_unavailable')
                self.assertNotIn('event_spike_counts', b)
                self.assertEqual(f['comparisons/0/u/head/onset'].attrs['status'], 'timing_unavailable')

    def test_overlap_rejected_for_independent_condition_comparison(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            inputs = self.create_inputs(root)
            with h5py.File(inputs['responses'], 'a') as f:
                del f['conditions/sound_2']
                f.copy('conditions/sound_1', 'conditions/sound_2')
            with self.assertRaisesRegex(ValueError, 'shares events'):
                run_session(inputs, root/'bad.h5', root/'bad.csv', root/'bad2.csv', {}, 'test')


if __name__ == '__main__':
    unittest.main()
