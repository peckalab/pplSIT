"""Regression tests for cache completeness, safe reuse, and stage boundaries."""
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow'))
from utils.sound_responses.settings import Settings, fitting_config, timing_config
from utils.sound_responses.segmentation import fit_spikes
from utils.sound_responses.bayesian_blocks import estimate_timing, select_timing
from utils.sound_responses.fits_session import fit_session
from utils.sound_responses.session import select_session


class CacheTests(unittest.TestCase):
    @staticmethod
    def inputs(root):
        events = np.column_stack((np.arange(40) * .25, np.ones(40)))
        parameters = {'sound': {'latency': .25, 'sounds': {
            'background': {'duration': .05, 'freq': 1000, 'channels': [3, 8]}}}}
        with h5py.File(root / 'meta.h5', 'w') as f:
            group = f.create_group('processed')
            group.attrs['parameters'] = json.dumps(parameters)
            group.create_dataset('timeline', data=[[0, 0], [10, 0]])
            group.create_dataset('sound_events', data=events)
        with h5py.File(root / 'units.h5', 'w') as f:
            rng = np.random.default_rng(42)
            f.create_dataset('unit-a/spike_times', data=np.sort(rng.uniform(0, 10, 200)))
            f.create_dataset('unit-empty/spike_times', data=[])

    def test_config_projections_isolate_timing_changes(self):
        base = {'window_s': .24, 'backend': 'gpu'}
        changes = {'offset_search_relative_to_duration_s': [-.007, .013],
                   'short_latency_s': .01, 'edge_margin_s': .002, 'min_block_s': .002,
                   'psth_bin_s': .01}
        for key, value in changes.items():
            changed = {**base, key: value}
            self.assertEqual(fitting_config(base), fitting_config(changed))
            self.assertNotEqual(timing_config(base), timing_config(changed))
        for key, value in {'window_s': .2, 'p0_values': [.05, .001],
                           'time_ranges_s': [[1, 9]], 'unit_names': ['unit-a'],
                           'conditions': {'bgr': [1]}, 'min_events': 30}.items():
            self.assertNotEqual(fitting_config(base), fitting_config({**base, key: value}))
            self.assertEqual(timing_config(base), timing_config({**base, key: value}))

    def test_cache_keeps_levels_beyond_original_early_stop(self):
        settings = Settings(p0_values=(.05, .01))
        spikes = np.linspace(.01, .23, 30)
        edges = {.05: np.array([-.25, .0005, .5]), .01: np.array([-.25, .02, .08, .5])}
        calls = []
        def backend(data, p0):
            calls.append(p0)
            return edges[p0]
        cached = fit_spikes(spikes, 30, settings, backend)
        self.assertEqual(calls, [.05, .01])
        original = select_timing(cached['fits'], 30, 30, .05, settings)
        self.assertEqual(original['is_modulated'], 0)
        self.assertEqual(len(original['fits']), 1)
        changed = replace(settings, edge_margin_s=.0001)
        selected = select_timing(cached['fits'], 30, 30, .05, changed)
        self.assertEqual(selected['onset_latency_s'], .02)
        # Offset search still stops at the first penalty without a candidate.
        self.assertEqual(selected['offset_status'], 'not_detected')
        self.assertEqual(len(selected['fits']), 2)
        reference = estimate_timing(spikes, 30, .05, changed, backend)
        for key in ('onset_latency_s', 'offset_latency_s', 'modulation_p0', 'offset_p0'):
            np.testing.assert_equal(selected[key], reference[key])

    def test_unused_later_numerical_failure_does_not_change_early_result(self):
        settings = Settings(p0_values=(.05, .01))
        def backend(data, p0):
            if p0 == .01:
                raise FloatingPointError('later level failed')
            return np.array([-.25, .5])
        cached = fit_spikes(np.linspace(.01, .23, 30), 30, settings, backend)
        self.assertEqual(cached['fits'][1]['status'], 'backend_numerical_failure')
        result = select_timing(cached['fits'], 30, 30, .05, settings)
        self.assertEqual((result['status'], result['is_modulated']), ('analyzed', 0))

    def test_reselection_and_rebinning_never_refit_or_modify_cache(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            self.inputs(root)
            config = {'window_s': .24, 'include_pooled': False}
            # Every penalty has a plausible onset/offset pair; count backend calls.
            with patch('hepstats.modeling.bayesian_blocks', return_value=np.array([-.24, .02, .075, .48])) as backend:
                fit_session(root / 'meta.h5', root / 'units.h5', root / 'fits.h5', config)
                self.assertEqual(backend.call_count, 5)
            before = hashlib.sha256((root / 'fits.h5').read_bytes()).hexdigest()
            modified = (root / 'fits.h5').stat().st_mtime_ns
            with patch('hepstats.modeling.bayesian_blocks', side_effect=AssertionError('must not refit')), \
                 patch('utils.sound_responses.gpu.GPUBackend', side_effect=AssertionError('must not use GPU')):
                rows = select_session(root / 'fits.h5', root / 'new.h5', root / 'new.csv', config)
                changed = select_session(root / 'fits.h5', root / 'old.h5', root / 'old.csv',
                                         {**config, 'offset_search_relative_to_duration_s': [-.007, .013], 'psth_bin_s': .01})
            self.assertEqual(rows[0]['onset_latency_s'], changed[0]['onset_latency_s'])
            self.assertEqual(rows[0]['offset_latency_s'], .075)
            with h5py.File(root / 'new.h5') as new, h5py.File(root / 'old.h5') as old:
                a, b = new['conditions/sound_1/units/unit-a'], old['conditions/sound_1/units/unit-a']
                self.assertNotEqual(len(a['psth_counts']), len(b['psth_counts']))
                self.assertEqual(sum(a['psth_counts']), sum(b['psth_counts']))
                self.assertEqual(new.attrs['selection_backend'], 'cached_edges_cpu')
                self.assertEqual(json.loads(new.attrs['segmentation_cache'])['sha256'], before)
            self.assertEqual(before, hashlib.sha256((root / 'fits.h5').read_bytes()).hexdigest())
            self.assertEqual(modified, (root / 'fits.h5').stat().st_mtime_ns)

    def test_rejects_wrong_fit_settings_incomplete_cache_and_inplace_overwrite(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            self.inputs(root)
            config = {'window_s': .24, 'include_pooled': False}
            with patch('hepstats.modeling.bayesian_blocks', return_value=np.array([-.24, .48])):
                fit_session(root / 'meta.h5', root / 'units.h5', root / 'fits.h5', config)
            for key, value in [('window_s', .2), ('p0_values', [.05]), ('min_spikes', 30),
                               ('time_ranges_s', [[1, 9]]), ('backend', 'gpu')]:
                with self.assertRaisesRegex(ValueError, 'regenerate fits.h5'):
                    select_session(root / 'fits.h5', root / 'out.h5', root / 'out.csv', {key: value})
            with self.assertRaisesRegex(ValueError, 'paths must be different'):
                select_session(root / 'fits.h5', root / 'fits.h5', root / 'out.csv', {})
            with h5py.File(root / 'fits.h5', 'r+') as f:
                del f['conditions/sound_1/units/unit-a/fits/4']
            with self.assertRaisesRegex(ValueError, 'incomplete penalty cache'):
                select_session(root / 'fits.h5', root / 'out.h5', root / 'out.csv', {})
            with h5py.File(root / 'fits.h5', 'r+') as f:
                f.attrs['complete'] = False
            with self.assertRaisesRegex(ValueError, 'complete fits.h5'):
                select_session(root / 'fits.h5', root / 'out.h5', root / 'out.csv', {})
            with h5py.File(root / 'legacy.h5', 'w') as f:
                f.attrs['schema_version'] = 1
            with self.assertRaisesRegex(ValueError, 'legacy responses.h5'):
                select_session(root / 'legacy.h5', root / 'out.h5', root / 'out.csv', {})


if __name__ == '__main__':
    unittest.main()
