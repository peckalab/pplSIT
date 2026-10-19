"""Run with the sound_responses_bayesian_blocks environment; no recording files needed."""
import json
from dataclasses import replace
from pathlib import Path
import sys
import tempfile
import unittest

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow'))
from utils.sound_responses.bayesian_blocks import (Settings, align_spikes, diagnostic_psth,
                                                  estimate_timing, select_events)
from utils.sound_responses.session import make_conditions, run_session
from utils.sound_responses.plots import plot_results


class EventTests(unittest.TestCase):
    def test_complete_windows_and_exclusion_reasons(self):
        events = np.column_stack(([0, .25, .49, .75, 1.0], [1, 2, 1, 2, 1]))
        valid, flags = select_events(events, .25, (0, 1), [[0, .5], [.75, 1]])
        np.testing.assert_array_equal(valid, [True, False, False, True, False])
        np.testing.assert_array_equal(flags, [0, 2, 4, 0, 5])

    def test_alignment_is_half_open_and_sorts_spikes(self):
        a = align_spikes([.5, .25, .249, 0, -.1, .6], [0, .25], .25)
        np.testing.assert_allclose(a, [0, .249, 0])
        edges, counts, rate = diagnostic_psth(a, 2, Settings())
        self.assertEqual(sum(counts), 3)
        self.assertAlmostEqual(np.sum(rate * np.diff(edges)), 1.5)

    def test_bad_event_inputs_fail(self):
        for events in ([[0, 1], [0, 2]], [[1, 1], [0, 2]], [[np.nan, 1]], [[0, 1.5]]):
            with self.assertRaises(ValueError):
                select_events(events, .25, (0, 2))
        with self.assertRaises(ValueError):
            select_events([[0, 1]], .25, (0, 2), [[0, 1], [.5, 2]])

    def test_conditions_include_sync_channel_metadata_without_treating_it_as_sound(self):
        events = np.array([[0, 1], [.25, 2], [.5, 0], [.75, -1]])
        catalog = {'1': {'name': 'background', 'duration': .05, 'channels': [3, 8]},
                   '2': {'name': 'target', 'duration': .05, 'channels': [3, 8]}}
        conditions = make_conditions(events, catalog)
        self.assertEqual([c['sound_ids'] for c in conditions], [[1], [2], [1, 2]])
        self.assertEqual(conditions[-1]['duration_s'], .05)
        catalog['2']['duration'] = .1
        pooled = make_conditions(events, catalog)[-1]
        self.assertTrue(np.isnan(pooled['duration_s']))
        self.assertEqual(pooled['duration_status'], 'mixed_durations')
        with self.assertRaises(ValueError):
            make_conditions(events, catalog, {'absent': [9]})


class TimingTests(unittest.TestCase):
    @staticmethod
    def backend(sequence):
        iterator = iter(sequence)
        return lambda data, p0: np.array(next(iterator))

    def test_previous_successful_level_and_duration_dependent_offset(self):
        backend = self.backend([[-.25, .020, .078, .50], [-.25, .025, .079, .50], [-.25, .50]])
        r = estimate_timing(np.linspace(.001, .24, 30), 30, .050, Settings(), backend)
        self.assertEqual(r['is_modulated'], 1)
        self.assertEqual(r['onset_latency_s'], .020)
        self.assertEqual(r['modulation_p0'], .01)
        self.assertEqual(r['onset_edge_p0'], .05)
        self.assertEqual(r['offset_latency_s'], .079)
        self.assertEqual(r['offset_p0'], .01)
        # The same boundary can qualify for 50 ms sounds but not for 30 ms sounds.
        b = self.backend([[-.25, .020, .078, .50], [-.25, .50]])
        r = estimate_timing(np.linspace(.001, .24, 30), 30, .050, Settings(), b)
        self.assertEqual(r['offset_latency_s'], .078)
        b = self.backend([[-.25, .020, .078, .50], [-.25, .50]])
        r = estimate_timing(np.linspace(.001, .24, 30), 30, .030, Settings(), b)
        self.assertEqual(r['offset_status'], 'not_detected')

    def test_offset_search_shift_preserves_onset_and_uses_strict_bounds(self):
        spikes = np.linspace(.001, .24, 30)
        settings = Settings()
        old = replace(settings, offset_search_relative_to_duration_s=(-.007, .013))
        # Exact binary fractions keep endpoint checks independent of rounding.
        onset, duration = .03125, .0625
        low, high = onset + duration + np.asarray(settings.offset_search_relative_to_duration_s)
        for edge, old_detected, new_detected in [
                (low - .004, True, False), (low, True, False),
                (low + .005, True, True), (low + .017, False, True),
                (high, False, False)]:
            with self.subTest(edge=edge):
                fits = [[-.25, onset, edge, .5], [-.25, .5]]
                before = estimate_timing(spikes, 30, duration, old, self.backend(fits))
                after = estimate_timing(spikes, 30, duration, settings, self.backend(fits))
                for key in ('onset_latency_s', 'onset_edge_p0', 'modulation_p0', 'is_modulated'):
                    self.assertEqual(before[key], after[key])
                self.assertEqual(before['offset_status'] == 'detected', old_detected)
                self.assertEqual(after['offset_status'] == 'detected', new_detected)

    def test_offset_search_must_fit_inside_observation_window(self):
        fits = [[-.25, .17, .225, .5], [-.25, .5]]
        r = estimate_timing(np.linspace(.001, .23, 30), 30, .05,
                            Settings(window_s=.24), self.backend(fits))
        self.assertEqual(r['onset_status'], 'detected')
        self.assertEqual(r['offset_status'], 'search_outside_window')
        self.assertTrue(np.isnan(r['offset_latency_s']))

    def test_dated_short_block_filter_removes_right_edge(self):
        b = self.backend([[-.25, .02, .0205, .08, .5], [-.25, .5]])
        r = estimate_timing(np.linspace(.001, .24, 30), 30, .057, Settings(), b)
        np.testing.assert_array_equal(r['fits'][0]['edges_s'], [-.25, .02, .08, .5])
        self.assertEqual(r['onset_latency_s'], .02)

    def test_short_latency_adjustment_and_outside_window_guard(self):
        for edges, expected in [([-.25, .004, .020, .080, .5], .020),
                                ([-.25, .004, .5], np.nan)]:
            b = self.backend([edges, [-.25, .5]])
            r = estimate_timing(np.linspace(.001, .24, 30), 30, .057, Settings(), b)
            self.assertTrue(r['short_latency_adjusted'])
            self.assertEqual(r['is_modulated'], 1)
            if np.isnan(expected):
                self.assertTrue(np.isnan(r['onset_latency_s']))
                self.assertEqual(r['onset_status'], 'reference_boundary_outside_window')
            else:
                self.assertEqual(r['onset_latency_s'], expected)

    def test_unavailable_is_not_negative_modulation(self):
        for spikes, n, settings, status in [([], 30, Settings(), 'insufficient_spikes'),
                                           ([.01] * 30, 2, Settings(), 'insufficient_events'),
                                           (np.linspace(.01, .2, 40), 30, Settings(max_spikes=30), 'spike_limit_exceeded')]:
            r = estimate_timing(spikes, n, .057, settings)
            self.assertEqual((r['status'], r['is_modulated']), (status, -1))
        b = self.backend([[-.25, .5]])
        r = estimate_timing(np.linspace(.001, .24, 30), 30, .057, Settings(), b)
        self.assertEqual(r['is_modulated'], 0)
        self.assertEqual(r['status'], 'analyzed')

    def test_real_backend_flat_rate_has_no_internal_boundary(self):
        spikes = (np.arange(400) + .5) * .25 / 400
        r = estimate_timing(spikes, 80, .057, Settings())
        self.assertEqual(r['is_modulated'], 0)

    def test_backend_numerical_failure_is_explicitly_unavailable(self):
        # Subtraction produces almost-coincident phases for a periodic spike train.
        spikes = align_spikes((np.arange(2000) + .5) / 100, np.arange(80) * .25, .25)
        r = estimate_timing(spikes, 80, .057, Settings())
        self.assertEqual(r['status'], 'backend_numerical_failure')
        self.assertEqual(r['is_modulated'], -1)
        self.assertTrue(r['diagnostic'])

    def test_real_backend_detects_known_evoked_interval(self):
        rng = np.random.default_rng(712)
        spikes = np.r_[rng.uniform(0, .25, 160), rng.uniform(.035, .092, 1000)]
        r = estimate_timing(spikes, 80, .050, Settings())
        self.assertEqual(r['is_modulated'], 1)
        self.assertAlmostEqual(r['onset_latency_s'], .035, delta=.004)
        self.assertAlmostEqual(r['offset_latency_s'], .092, delta=.004)

    def test_frozen_legacy_backend_and_dated_reference_results(self):
        from hepstats.modeling import bayesian_blocks
        fixture = json.loads((Path(__file__).parent/'fixtures/sound_responses_bayesian_blocks_reference.json').read_text())
        for case in fixture['cases']:
            with self.subTest(kind=case['kind']):
                rng = np.random.default_rng(case['seed'])
                if case['kind'] == 'flat':
                    a = (np.arange(400) + .5) * .25 / 400
                elif case['kind'] == 'excited':
                    a = np.r_[rng.uniform(0, .25, 160), rng.uniform(.035, .092, 1000)]
                else:
                    a = rng.uniform(0, .25, 1600)
                    a = a[(a < .035) | (a > .092)]
                wrapped = np.r_[a - .25, a, a + .25]
                for p0, expected in zip(Settings().p0_values, case['edges']):
                    np.testing.assert_allclose(bayesian_blocks(wrapped, p0=p0), expected, rtol=0, atol=1e-12)
                # Original offset code hardcodes onset + 50–70 ms for 50 ms sounds.
                result = estimate_timing(a, 80, .050, Settings())
                self.assertEqual(result['is_modulated'], int(case['is_modulated']))
                for key in ('onset_latency_s', 'modulation_p0', 'offset_latency_s', 'offset_p0'):
                    expected = case[key] if case[key] is not None else np.nan
                    np.testing.assert_allclose(result[key], expected, rtol=0, atol=1e-12, equal_nan=True)

    def test_invalid_settings(self):
        for kwargs in ({'window_s': 0}, {'p0_values': [.01, .05]}, {'p0_values': [np.nan]},
                       {'min_events': 1.5}, {'max_spikes': 2}, {'offset_search_relative_to_duration_s': [0, 0]}):
            with self.assertRaises(ValueError):
                Settings(**kwargs)
        with self.assertRaises(ValueError):
            Settings.from_config({}, {})


class SessionTests(unittest.TestCase):
    def test_hdf_csv_pdf_and_passive_metadata_fallback(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            events = np.column_stack((np.arange(80) * .25, np.tile([1, 2], 40)))
            parameters = {'sound': {'latency': .25, 'sounds': {
                'F_1000': {'duration': .05, 'freq': 1000, 'channels': [3, 8]},
                'F_2000': {'duration': .1, 'freq': 2000, 'channels': [3, 8]}}}}
            with h5py.File(root/'meta.h5', 'w') as f:
                g = f.create_group('processed'); g.attrs['parameters'] = json.dumps(parameters)
                g.create_dataset('timeline', data=[[0, 0], [20, 0]])
                g.create_dataset('sound_events', data=events)
            with h5py.File(root/'units.h5', 'w') as f:
                f.create_dataset('unit-a/spike_times', data=(np.arange(2000) + .5) / 100)
                f.create_dataset('unit-empty/spike_times', data=[])
            rows = run_session(root/'meta.h5', root/'units.h5', root/'responses.h5', root/'summary.csv', {})
            self.assertEqual(len(rows), 6)
            with h5py.File(root/'responses.h5', 'r') as f:
                self.assertEqual(f.attrs['method'], 'amaro_07062021_bayesian_blocks')
                self.assertIn('not a responsiveness-test p-value', f.attrs['p0_interpretation'])
                self.assertEqual(f['conditions/sound_1'].attrs['label'], 'F_1000')
                self.assertEqual(f['conditions/pooled'].attrs['duration_status'], 'mixed_durations')
                self.assertEqual(f['conditions/sound_1/units/unit-empty'].attrs['is_modulated'], -1)
                self.assertEqual(len(f['conditions/sound_1/event_indices']), 40)
            plot_results(root/'responses.h5', root/'figures/responses.pdf')
            self.assertTrue((root/'figures/responses.pdf').read_bytes().startswith(b'%PDF'))
            self.assertEqual(len((root/'summary.csv').read_text().splitlines()), 7)

    def test_empty_conditions_write_explicit_report(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            with h5py.File(root/'meta.h5', 'w') as f:
                g = f.create_group('processed');g.attrs['parameters'] = json.dumps({'sound': {'latency': .25}})
                g.create_dataset('timeline', data=[[0, 0], [1, 0]])
                g.create_dataset('sound_events', data=[[0, 0], [.25, 0]])
            with h5py.File(root/'units.h5', 'w') as f:
                f.create_dataset('empty/spike_times', data=[])
            rows = run_session(root/'meta.h5', root/'units.h5', root/'r.h5', root/'r.csv', {})
            self.assertEqual(rows, [])
            plot_results(root/'r.h5', root/'r.pdf')
            self.assertEqual(len((root/'r.csv').read_text().splitlines()), 1)


if __name__ == '__main__':
    unittest.main()
