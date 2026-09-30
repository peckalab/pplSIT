"""Timestamp repair regression tests; no recording files or Snakemake needed."""
import os
import ast
from pathlib import Path
import sys
import tempfile
import unittest
import warnings

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow'))
from utils import sync


class TimestampTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def save(self, values, folder='stream', samples=True):
        folder = self.root / folder
        folder.mkdir(parents=True, exist_ok=True)
        path = folder / 'timestamps.npy'
        np.save(path, values)
        if samples:
            np.save(folder / 'sample_numbers.npy', np.arange(len(values), dtype=np.int64))
        return str(path)

    def test_clean_unchanged(self):
        a = 10 + np.arange(20001) / 30000
        path = self.save(a)
        out, _ = sync._load_timestamps_for_sync(path, 'test')
        np.testing.assert_array_equal(out, a)

    def test_invalid_values_rejected(self):
        for a in ([0, 1, np.nan, 3], [0, 1, np.inf], [-1, -1, -1], [], [1]):
            with self.subTest(values=a), self.assertRaises(ValueError):
                sync._load_timestamps_for_sync(self.save(a), 'test')

    def test_single_or_multiple_leading_invalid_values(self):
        for count in (1, 4):
            a = 10 + np.arange(100) / 30000
            broken = a.copy()
            broken[:count] = -1
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                out, _ = sync._load_timestamps_for_sync(self.save(broken), 'test')
            np.testing.assert_allclose(out, a, rtol=0, atol=1e-12)
            np.testing.assert_array_equal(out[count:], broken[count:])

    def test_midstream_negative_rejected(self):
        a = 10 + np.arange(20001) / 30000
        a[10000] = -1
        with self.assertRaises(ValueError):
            sync._load_timestamps_for_sync(self.save(a), 'test')

    def test_bounded_inversion_preserves_surroundings(self):
        a = 10 + np.arange(20001) / 30000
        a[10000] = a[9999] - 15 / 30000
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            out, _ = sync._load_timestamps_for_sync(self.save(a), 'test')
        self.assertTrue(np.all(np.diff(out) > 0))
        np.testing.assert_array_equal(out[:10000], a[:10000])
        np.testing.assert_array_equal(out[10001:], a[10001:])

    def test_leading_repair_does_not_hide_large_jump(self):
        a = 10 + np.arange(30001) / 30000
        a[:4] = -1
        a[15000:] -= 0.333
        with warnings.catch_warnings(), self.assertRaises(ValueError):
            warnings.simplefilter('ignore')
            sync._load_timestamps_for_sync(self.save(a), 'test')

    def test_noncontiguous_samples_rejected(self):
        a = 10 + np.arange(20001) / 30000
        a[10000] = a[9999] - 15 / 30000
        path = self.save(a)
        samples = np.arange(len(a)); samples[10000:] += 1
        np.save(Path(path).with_name('sample_numbers.npy'), samples)
        with self.assertRaises(ValueError):
            sync._load_timestamps_for_sync(path, 'test')

    def test_default_keeps_legacy_inversion_result(self):
        a = 10 + np.arange(20001) / 30000
        a[10000] = a[9999] - 0.5 / 30000
        path = self.save(a)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            expected = sync._repair_legacy_timestamp_inversions(a, path, 'ADC')
            out, _ = sync._load_adc_timestamps_for_sync(path, None, 30000)
        np.testing.assert_array_equal(out, expected)
        np.testing.assert_array_equal(sync._load_ephys_timestamps_for_sync(path), a)

    def test_default_also_rejects_nonfinite(self):
        for a in ([0, 1, np.nan, 3], [0, 1, np.inf]):
            with self.subTest(a=a), self.assertRaises(ValueError):
                sync._load_adc_timestamps_for_sync(self.save(a), None, 30000)

    def test_shared_origin_with_units_and_subset_of_probes(self):
        # Execute the actual units helper without running the Snakemake script.
        path = Path(__file__).resolve().parents[1] / 'workflow/scripts/units.py'
        node = next(n for n in ast.parse(path.read_text()).body
                    if isinstance(n, ast.FunctionDef) and n.name == 'get_session_ephys_t0')
        env = dict(os=os, np=np, reference_ephys_timestamp_path=sync.reference_ephys_timestamp_path,
                   _load_ephys_timestamps_for_sync=sync._load_ephys_timestamps_for_sync)
        exec(compile(ast.Module(body=[node], type_ignores=[]), str(path), 'exec'), env)
        a = 10 + np.arange(100) / 30000
        broken = a.copy(); broken[:4] = -1
        self.save(broken, 'ephys/ProbeA')
        self.save(a + 1, 'ephys/ProbeB')
        root = str(self.root / 'ephys')
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            unit_t0 = env['get_session_ephys_t0'](root, ['ProbeB'], True)
            sound_t0 = sync._load_ephys_timestamps_for_sync(sync.reference_ephys_timestamp_path(root), True)[0]
        self.assertEqual(unit_t0, sound_t0)
        self.assertAlmostEqual(unit_t0, 10)
        self.assertEqual(env['get_session_ephys_t0'](root, ['ProbeA', 'ProbeB']), -1)

    def test_original_sample_discovery_does_not_write(self):
        a = 10 + np.arange(100) / 30000; a[:4] = -1
        original = self.save(a, 'session/original/ProbeA')
        staged = self.root / 'session/ephys/ProbeA'
        staged.mkdir(parents=True)
        os.link(original, staged / 'timestamps.npy')
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            out, _ = sync._load_timestamps_for_sync(str(staged / 'timestamps.npy'), 'test')
        self.assertAlmostEqual(out[0], 10)
        self.assertFalse((staged / 'sample_numbers.npy').exists())
        np.testing.assert_array_equal(np.load(original), a)

    def test_missing_samples_blocks_leading_repair(self):
        with self.assertRaises(FileNotFoundError):
            sync._load_timestamps_for_sync(self.save([-1, 10, 11], samples=False), 'test')

    def test_unbounded_and_unanchored_inversions_rejected(self):
        for kind in ('large', 'end', 'dense', 'span'):
            a = 10 + np.arange(30001) / 30000
            if kind == 'large': a[15000] = a[14999] - 100 / 30000
            if kind == 'end': a[-1] = a[-2]
            if kind == 'dense': a[100:120] = a[99]
            if kind == 'span': a[15000:15070] = a[14999] - 1 / 30000 + np.arange(70) * 1e-9
            with self.subTest(kind=kind), self.assertRaises(ValueError):
                sync._load_timestamps_for_sync(self.save(a), 'test')

    def test_clean_sound_sync_outputs_unchanged(self):
        fs = 1000
        times = 10 + np.arange(4000) / fs
        path = self.save(times)
        adc = np.zeros(4000, dtype=np.int16)
        for index in (500, 1500, 2500):
            adc[index:index + 50] = 1000
        adc.tofile(self.root / 'adc.dat')
        np.savetxt(self.root / 'events.csv', [[100, 0], [104, 1]], delimiter=',', header='time,event')
        np.savetxt(self.root / 'sounds.csv', [[100.5, 1], [101.5, 2], [102.5, 1]], delimiter=',', header='time,sound')
        kwargs = dict(adc_file=str(self.root / 'adc.dat'), adc_ts_file=path,
                      ephys_ts_file=path, events_file=str(self.root / 'events.csv'),
                      sounds_file=str(self.root / 'sounds.csv'), channel=0, ch_no=1)
        for sample_rate in (fs, None):
            legacy = sync.get_sound_events_from_ADC(**kwargs, s_rate=sample_rate)
            repaired = sync.get_sound_events_from_ADC(**kwargs, s_rate=sample_rate, repair_timestamps=True)
            self.assertEqual(len(legacy[0]), 3)
            for expected, actual in zip(legacy, repaired):
                np.testing.assert_array_equal(actual, expected)

    def test_sample_count_mismatch_rejected(self):
        path = self.save([-1, -1, 10, 11])
        np.save(Path(path).with_name('sample_numbers.npy'), np.arange(3))
        with self.assertRaises(ValueError):
            sync._load_timestamps_for_sync(path, 'test')


if __name__ == '__main__':
    unittest.main()
