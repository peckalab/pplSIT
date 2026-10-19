"""Shared settings and stage-specific configuration projections."""
from dataclasses import dataclass, fields, asdict
import numpy as np

@dataclass(frozen=True)
class Settings:
    window_s: float = 0.25
    p0_values: tuple = (0.05, 0.01, 0.001, 0.0001, 0.00001)
    min_block_s: float = 0.001
    edge_margin_s: float = 0.001
    short_latency_s: float = 0.007
    offset_search_relative_to_duration_s: tuple = (0.0, 0.020)
    min_events: int = 20
    min_spikes: int = 20
    max_spikes: int | None = None
    psth_bin_s: float = 0.005

    def __post_init__(self):
        for name in ('window_s', 'min_block_s', 'edge_margin_s', 'short_latency_s', 'psth_bin_s'):
            value = getattr(self, name)
            if not np.isfinite(value) or value <= 0:
                raise ValueError(f'{name} must be finite and positive')
        if self.window_s <= 2 * self.edge_margin_s or self.short_latency_s >= self.window_s:
            raise ValueError('Analysis window is too short for the configured boundary criteria')
        if self.min_block_s >= self.window_s:
            raise ValueError('min_block_s must be smaller than window_s')
        p0 = np.asarray(self.p0_values, dtype=float)
        if (p0.ndim != 1 or not len(p0) or not np.all(np.isfinite(p0))
                or np.any((p0 <= 0) | (p0 >= 1)) or np.any(np.diff(p0) >= 0)):
            raise ValueError('p0_values must be strictly decreasing values between zero and one')
        offsets = np.asarray(self.offset_search_relative_to_duration_s, dtype=float)
        if offsets.shape != (2,) or not np.all(np.isfinite(offsets)) or offsets[0] >= offsets[1]:
            raise ValueError('offset_search_relative_to_duration_s must be an increasing pair')
        for name in ('min_events', 'min_spikes', 'max_spikes'):
            value = getattr(self, name)
            if name == 'max_spikes' and value is None:
                continue
            if isinstance(value, bool) or not isinstance(value, (int, np.integer)) or value < 1:
                raise ValueError(f'{name} must be a positive integer')
        if self.max_spikes is not None and self.max_spikes < self.min_spikes:
            raise ValueError('max_spikes must be >= min_spikes')

    @classmethod
    def from_config(cls, config, parameters):
        values = {f.name: config[f.name] for f in fields(cls) if f.name in config}
        if config.get('window_s') is None:
            window = parameters.get('sound', {}).get('latency')
            if window is None:
                raise ValueError('Set window_s explicitly: session sound.latency is missing')
            values['window_s'] = float(window)
        return cls(**values)



FIT_SETTING_KEYS = ('window_s', 'p0_values', 'min_events', 'min_spikes', 'max_spikes')
TIMING_SETTING_KEYS = ('min_block_s', 'edge_margin_s', 'short_latency_s',
                      'offset_search_relative_to_duration_s', 'psth_bin_s')


def fitting_config(config):
    """Only values that affect fitted data, QC, segmentation, or execution."""
    defaults = asdict(Settings())
    result = {key: config.get(key, None if key == 'window_s' else defaults[key])
              for key in FIT_SETTING_KEYS}
    result.update(backend=config.get('backend', 'cpu'), conditions=config.get('conditions'),
                  include_pooled=config.get('include_pooled', True),
                  unit_names=config.get('unit_names'), time_ranges_s=config.get('time_ranges_s'))
    if result['backend'] == 'gpu':
        result['gpu_memory_mb'] = config.get('gpu_memory_mb', 4096)
    return result


def timing_config(config):
    defaults = asdict(Settings())
    return {key: config.get(key, defaults[key]) for key in TIMING_SETTING_KEYS}
