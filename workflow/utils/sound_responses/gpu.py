"""Single-GPU, unweighted float64 Bayesian blocks matching hepstats 0.9.2.

CUDA changes only segmentation. The shared estimate_timing function still
selects boundaries, applies QC, and records only CPU-equivalent penalty levels.
CuPy is imported lazily so CPU-only installations do not need CUDA.
"""
from contextlib import contextmanager
import fcntl
import os
from pathlib import Path
import time

import numpy as np

from .segmentation import prepare_input, fit_spikes, empty_segmentation


def _positive_integer(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)) or value < 1:
        raise ValueError(f'{name} must be a positive integer')
    return int(value)


def _prepare(data):
    data, counts = np.unique(np.asarray(data, dtype=np.float64), return_counts=True)
    if len(data) < 2 or not np.all(np.isfinite(data)):
        raise ValueError('GPU fitting requires at least two distinct, finite timestamps')
    if len(data) > np.iinfo(np.int32).max or counts.sum() >= 2**53:
        raise ValueError('GPU fit exceeds exact index/count representation')
    edges = np.r_[data[:1], .5 * (data[1:] + data[:-1]), data[-1:]]
    lengths = data[-1] - edges
    if np.any(np.diff(lengths) >= 0) or not np.all(np.isfinite(lengths)):
        raise FloatingPointError('nonpositive Bayesian-block cell width')
    return edges, lengths, np.r_[0., np.cumsum(counts.astype(np.float64))]


@contextmanager
def _device_lock(cp):
    # Stable host-local path even when Snakemake gives each job a different TMPDIR.
    # Also serializes this user's jobs from separate Snakemake invocations.
    bus = cp.cuda.Device(0).pci_bus_id.replace(':', '_').replace('.', '_')
    path = Path('/tmp') / f'pplsit-sound-responses-gpu-{os.getuid()}-{bus}.lock'
    with path.open('a') as lock:
        print(f'Waiting for Bayesian-block GPU lock: {bus}', flush=True)
        fcntl.flock(lock, fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(lock, fcntl.LOCK_UN)


class GPUBackend:
    def __init__(self, settings, memory_mb=4096):
        self.settings = settings
        self.limit = _positive_integer(memory_mb, 'gpu_memory_mb') * 1024**2
        if len(settings.p0_values) > 16:
            raise ValueError('GPU backend supports at most 16 p0 values; use CPU for longer sequences')
        visible = os.environ.get('CUDA_VISIBLE_DEVICES', '')
        if len(visible.split(',')) != 1 or not visible.strip() or visible.strip() == '-1':
            raise ValueError('GPU backend requires exactly one device in CUDA_VISIBLE_DEVICES '
                             '(prefer its GPU UUID); no automatic GPU selection or CPU fallback')
        import cupy as cp
        if cp.cuda.runtime.getDeviceCount() != 1:
            raise RuntimeError('GPU backend requires exactly one visible CUDA device')
        self.cp = cp
        self.pool = cp.cuda.MemoryPool()
        self.pool.set_limit(size=self.limit)
        self.kernel = cp.RawKernel(Path(__file__).with_name('gpu_kernel.cu').read_text(),
                                   'bayesian_blocks_batch', options=(
                                       '--std=c++11', '--fmad=false',
                                       f'-DN_PENALTIES={len(settings.p0_values)}'))
        props = cp.cuda.runtime.getDeviceProperties(0)
        name = props['name']
        self.provenance = dict(name='cupy_float64', cupy_version=cp.__version__,
                               device_name=name.decode() if isinstance(name, bytes) else str(name),
                               visible_devices=visible, pci_bus_id=cp.cuda.Device(0).pci_bus_id,
                               runtime_version=cp.cuda.runtime.runtimeGetVersion(),
                               driver_version=cp.cuda.runtime.driverGetVersion(),
                               memory_limit_bytes=self.limit, threads=512,
                               compiler_options=list(self.kernel.options), batches=0,
                               fitted_cases=0, peak_pool_bytes=0)

    def _fit(self, prepared):
        cp = self.cp
        sizes = np.array([len(x[0]) - 1 for x in prepared], dtype=np.int32)
        offsets = np.r_[0, np.cumsum(sizes, dtype=np.int64)[:-1]]
        total = int(sizes.sum())
        # Scalar Python power/log order preserves the reference's prior arithmetic.
        p0s = self.settings.p0_values
        table = np.array([[4 - np.log(73.53 * p * ((r + 1)**-.478))
                           for r in range(int(max(sizes)))] for p in p0s])
        arrays = [cp.asarray(np.concatenate([x[1] for x in prepared])),
                  cp.asarray(np.concatenate([x[2] for x in prepared])),
                  cp.asarray(np.concatenate([table[:, :n] for n in sizes], axis=1)),
                  cp.asarray(offsets), cp.asarray(sizes)]
        best = cp.empty((len(p0s), total), dtype=cp.float64)
        last = cp.empty((len(p0s), total), dtype=cp.int32)
        invalid = cp.zeros(len(prepared), dtype=cp.int32)
        self.kernel((len(prepared),), (512,), (*arrays, np.int64(total), best, last, invalid))
        cp.cuda.Stream.null.synchronize()
        invalid, last = cp.asnumpy(invalid), cp.asnumpy(last)
        self.provenance['peak_pool_bytes'] = max(self.provenance['peak_pool_bytes'], self.pool.total_bytes())
        results = []
        for case, (base, n, item) in enumerate(zip(offsets, sizes, prepared)):
            if invalid[case]:
                results.append(FloatingPointError('nonpositive Bayesian-block cell width'))
                continue
            fits = {}
            for k, p0 in enumerate(p0s):
                index = int(n)
                indices = [index]
                while index:
                    previous = int(last[k, base + index - 1])
                    if not 0 <= previous < index:
                        raise RuntimeError('Invalid GPU Bayesian-block traceback')
                    index = previous
                    indices.append(index)
                fits[p0] = item[0][indices[::-1]]
            results.append(fits)
        self.provenance['batches'] += 1
        self.provenance['fitted_cases'] += len(prepared)
        return results

    def fit(self, cases):
        """Return all raw fits for (aligned_spikes, n_events) tuples.

        Batches obey a conservative device-memory budget. A single oversized
        case fails explicitly; it is never silently subsampled or sent to CPU.
        """
        started = time.perf_counter()
        results = [None] * len(cases)
        with _device_lock(self.cp), self.cp.cuda.Device(0), self.cp.cuda.using_allocator(self.pool.malloc):
            prepared, indices, cost = [], [], 0

            def finish():
                nonlocal prepared, indices, cost
                if not prepared:
                    return
                print(f'GPU Bayesian blocks: fitting batch of {len(prepared)} unit/condition cases', flush=True)
                fits = self._fit(prepared)
                for index, cached in zip(indices, fits):
                    def replay(data, p0):
                        if isinstance(cached, Exception):
                            raise cached
                        return cached[p0]
                    results[index] = fit_spikes(*cases[index], self.settings, backend=replay)
                prepared, indices, cost = [], [], 0
                self.pool.free_all_blocks()

            for index, (spikes, n_events) in enumerate(cases):
                a, unavailable = prepare_input(spikes, n_events, self.settings)
                if unavailable is not None:
                    results[index] = unavailable
                    continue
                try:
                    item = _prepare(np.r_[a - self.settings.window_s, a, a + self.settings.window_s])
                except FloatingPointError as exc:
                    results[index] = empty_segmentation('backend_numerical_failure', n_events, len(a), str(exc))
                    continue
                size = (len(item[0]) - 1) * (16 + 20 * len(self.settings.p0_values)) + 64
                # Leave room for allocator rounding/fragmentation; CUDA context is additional.
                if size > self.limit * .7:
                    raise MemoryError('One GPU fit exceeds gpu_memory_mb; increase the budget or select CPU')
                if cost + size > self.limit * .7:
                    finish()
                prepared.append(item)
                indices.append(index)
                cost += size
            finish()
        self.provenance['elapsed_seconds'] = time.perf_counter() - started
        return results

    def estimate(self, cases):
        """Convenience API; production rules save fit() results before selection."""
        from .bayesian_blocks import select_timing
        raw = self.fit([(a, n) for a, n, duration in cases])
        return [select_timing(result['fits'], result['n_events'], result['n_spikes'],
                              case[2], self.settings, status=result['status'],
                              diagnostic=result['diagnostic'])
                for result, case in zip(raw, cases)]
