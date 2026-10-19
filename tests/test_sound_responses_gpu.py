"""GPU tests are opt-in: RUN_PPLSIT_GPU_TESTS=1 CUDA_VISIBLE_DEVICES=<one GPU>."""
from dataclasses import replace
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow'))
from utils.sound_responses.bayesian_blocks import Settings, align_spikes, estimate_timing
from utils.sound_responses.gpu import GPUBackend, _prepare
from utils.sound_responses.session import run_session


def assert_results_equal(test, actual, expected):
    for key in expected:
        if key == 'fits':
            test.assertEqual(len(actual[key]), len(expected[key]))
            for a, b in zip(actual[key], expected[key]):
                test.assertEqual(a['p0'], b['p0'])
                for name in ('raw_edges_s', 'edges_s'):
                    np.testing.assert_array_equal(a[name], b[name])
        elif isinstance(expected[key], float):
            np.testing.assert_array_equal(actual[key], expected[key], err_msg=key)
        else:
            test.assertEqual(actual[key], expected[key], key)


class GPUConfigurationTests(unittest.TestCase):
    def test_explicit_single_device_required_before_importing_cupy(self):
        for value in ('', '0,1', '-1'):
            with self.subTest(visible=value), patch.dict(os.environ, CUDA_VISIBLE_DEVICES=value):
                with self.assertRaisesRegex(ValueError, 'exactly one device'):
                    GPUBackend(Settings())

    def test_invalid_budget_and_penalty_count(self):
        for value in (0, -1, 1.5, True):
            with self.assertRaises(ValueError):
                GPUBackend(Settings(), memory_mb=value)
        with self.assertRaisesRegex(ValueError, 'at most 16'):
            GPUBackend(Settings(p0_values=tuple(np.geomspace(.05, .00001, 17))))

    def test_cell_width_cancellation_rejected(self):
        with self.assertRaises(FloatingPointError):
            _prepare([0., np.nextafter(0., 1.), 1.])


@unittest.skipUnless(os.environ.get('RUN_PPLSIT_GPU_TESTS') == '1', 'requires an explicitly selected GPU')
class GPUEquivalenceTests(unittest.TestCase):
    def test_frozen_reference_all_raw_edges_and_timing(self):
        fixture = json.loads((Path(__file__).parent / 'fixtures/sound_responses_bayesian_blocks_reference.json').read_text())
        settings = Settings()
        cases, expected = [], []
        for case in fixture['cases']:
            rng = np.random.default_rng(case['seed'])
            if case['kind'] == 'flat':
                a = (np.arange(400) + .5) * .25 / 400
            elif case['kind'] == 'excited':
                a = np.r_[rng.uniform(0, .25, 160), rng.uniform(.035, .092, 1000)]
            else:
                a = rng.uniform(0, .25, 1600)
                a = a[(a < .035) | (a > .092)]
            cases.append((a, 80, .05))
            expected.append(estimate_timing(a, 80, .05, settings))
        gpu = GPUBackend(settings)
        for actual, reference in zip(gpu.estimate(cases), expected):
            assert_results_equal(self, actual, reference)

    def test_custom_penalties_repeats_and_memory_bounded_batches(self):
        rng = np.random.default_rng(19)
        spikes = np.r_[rng.uniform(0, .24, 400), [.03] * 40, [.08] * 80]
        for p0s in ((.01,), (.05, .001)):
            settings = Settings(window_s=.24, p0_values=p0s)
            case = (spikes, 40, .05)
            expected = estimate_timing(*case, settings)
            gpu = GPUBackend(settings, memory_mb=1)
            results = gpu.estimate([case] * 20)
            self.assertGreater(gpu.provenance['batches'], 1)
            self.assertLessEqual(gpu.provenance['peak_pool_bytes'], 1024**2)
            for actual in results:
                assert_results_equal(self, actual, expected)

    def test_qc_numerical_failures_and_oversized_input(self):
        settings = Settings()
        bad = align_spikes((np.arange(2000) + .5) / 100, np.arange(80) * .25, .25)
        rng = np.random.default_rng(20)
        good = rng.uniform(0, .25, 150)
        cases = [(bad, 80, .05), (good, 40, .05), ([], 30, .05), (good, 2, .05)]
        results = GPUBackend(settings).estimate(cases)
        self.assertEqual(results[0]['status'], 'backend_numerical_failure')
        self.assertTrue(results[0]['diagnostic'])
        for actual, case in zip(results[1:], cases[1:]):
            assert_results_equal(self, actual, estimate_timing(*case, settings))
        with self.assertRaisesRegex(MemoryError, 'One GPU fit'):
            GPUBackend(settings, memory_mb=1).estimate([(np.linspace(.001, .24, 20000), 40, .05)])
        limited = GPUBackend(replace(settings, max_spikes=30))
        self.assertEqual(limited.estimate([(good, 40, .05)])[0]['status'], 'spike_limit_exceeded')
        self.assertEqual(limited.provenance['fitted_cases'], 0)
        self.assertEqual(limited.estimate([]), [])

    def test_session_hdf5_csv_matches_cpu(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            events = np.column_stack((np.arange(60) * .25, np.tile([1, 2], 30)))
            parameters = {'sound': {'latency': .25, 'sounds': {
                'background': {'duration': .05, 'freq': 1000, 'channels': [3, 8]},
                'target': {'duration': .05, 'freq': 2000, 'channels': [3, 8]}}}}
            with h5py.File(root / 'meta.h5', 'w') as f:
                group = f.create_group('processed')
                group.attrs['parameters'] = json.dumps(parameters)
                group.create_dataset('timeline', data=[[0, 0], [15, 0]])
                group.create_dataset('sound_events', data=events)
            with h5py.File(root / 'units.h5', 'w') as f:
                rng = np.random.default_rng(21)
                f.create_dataset('unit-a/spike_times', data=np.sort(rng.uniform(0, 15, 500)))
                f.create_dataset('unit-empty/spike_times', data=[])
            for backend in ('cpu', 'gpu'):
                run_session(root / 'meta.h5', root / 'units.h5', root / f'{backend}.h5',
                            root / f'{backend}.csv', {'backend': backend, 'window_s': .24})
            self.assertEqual((root / 'cpu.csv').read_text(), (root / 'gpu.csv').read_text())
            with h5py.File(root / 'cpu.h5') as cpu, h5py.File(root / 'gpu.h5') as gpu:
                self.assertEqual(cpu.attrs['backend'], 'cpu')
                self.assertEqual(gpu.attrs['backend'], 'gpu')
                self.assertEqual(json.loads(gpu.attrs['backend_provenance'])['name'], 'cupy_float64')
                def compare(name, obj):
                    if isinstance(obj, h5py.Dataset):
                        np.testing.assert_array_equal(obj[:], cpu[name][:])
                    if name.startswith('conditions/'):
                        for key, value in obj.attrs.items():
                            if isinstance(value, str):
                                self.assertEqual(value, cpu[name].attrs[key])
                            else:
                                np.testing.assert_array_equal(value, cpu[name].attrs[key])
                gpu.visititems(compare)


if __name__ == '__main__':
    unittest.main()
