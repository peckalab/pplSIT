import json
import os
from pathlib import Path


def _runtime_dir(name):
    tmpdir = str(getattr(snakemake.resources, "tmpdir", os.environ.get("TMPDIR", "/tmp")))
    path = os.path.join(tmpdir, "bombcell-runtime", name)
    os.makedirs(path, exist_ok=True)
    return path


thread_count = max(1, int(getattr(snakemake, "threads", 1)))
os.environ["LOKY_MAX_CPU_COUNT"] = str(thread_count)
os.environ["JOBLIB_TEMP_FOLDER"] = _runtime_dir("joblib")
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["BLIS_NUM_THREADS"] = "1"
os.environ["NUMBA_NUM_THREADS"] = "1"
os.environ.setdefault("NUMBA_CACHE_DIR", _runtime_dir("numba"))
os.environ.setdefault("MPLCONFIGDIR", _runtime_dir("matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", _runtime_dir("xdg-cache"))
os.environ.setdefault("CACHECACHE_DIR", _runtime_dir("cachecache"))
os.environ.setdefault("BOMBCELL_CACHE_DIR", _runtime_dir("bombcell-cache"))

home = os.environ.get("HOME")
if not home or not os.access(home, os.W_OK):
    os.environ["HOME"] = _runtime_dir("home")

import bombcell


def _config_value(section, key, default=None):
    return snakemake.config.get(section, {}).get(key, default)


def _infer_n_channels(dat_path, base_n_channels, dtype_bytes=2):
    size = os.path.getsize(dat_path)
    candidates = []
    for value in (base_n_channels, base_n_channels + 1):
        if value not in candidates:
            candidates.append(value)

    for n_channels in candidates:
        if size % (dtype_bytes * n_channels) == 0:
            return n_channels

    raise ValueError(
        f"Cannot infer channel count for {dat_path}. "
        f"File size {size} is not divisible by int16 samples for {candidates} channels."
    )


def _required_outputs():
    return [
        str(snakemake.output.unit_type),
        str(snakemake.output.metrics_csv),
        str(snakemake.output.metrics_parquet),
    ]


ks_dir = os.path.dirname(str(snakemake.input.st))
raw_file = str(snakemake.input.dat)
save_path = os.path.dirname(str(snakemake.output.ready))
settings_path = str(snakemake.input.settings)

os.makedirs(save_path, exist_ok=True)

with open(settings_path, "r") as f:
    ks_settings = json.load(f)

kilosort_version = int(_config_value("bombcell", "kilosort_version", 4))
gain_to_uV = _config_value("bombcell", "gain_to_uV", None)

param = bombcell.get_default_parameters(
    ks_dir,
    raw_file=raw_file,
    kilosort_version=kilosort_version,
    gain_to_uV=gain_to_uV,
)

param_overrides = _config_value("bombcell", "param_overrides", {}) or {}
if not isinstance(param_overrides, dict):
    raise ValueError("bombcell.param_overrides must be a mapping of Bombcell parameter names to values.")
param.update(param_overrides)

base_n_channels = int(_config_value("bombcell", "n_channels", _config_value("lfp", "n_channels", 384)))
n_chan_bin = _infer_n_channels(raw_file, base_n_channels)

param["ephys_sample_rate"] = float(ks_settings["fs"])
param["nChannels"] = n_chan_bin
param["nSyncChannels"] = max(0, n_chan_bin - base_n_channels)
param["saveAsTSV"] = True
param["unit_type_for_phy"] = True
param["plotDetails"] = False
param["plotGlobal"] = False
param["splitGoodAndMua_NonSomatic"] = bool(
    _config_value("bombcell", "split_good_and_mua_non_somatic", True)
)
param["savePlots"] = bool(_config_value("bombcell", "save_plots", False))
param["extractRaw"] = bool(_config_value("bombcell", "extract_raw", True))
param["reextractRaw"] = bool(_config_value("bombcell", "reextract_raw", False))
param["decompress_data"] = bool(_config_value("bombcell", "decompress_data", False))
param["nRawSpikesToExtract"] = int(_config_value("bombcell", "n_raw_spikes_to_extract", 100))
param["verbose"] = bool(_config_value("bombcell", "verbose", True))

bombcell.run_bombcell(
    ks_dir,
    save_path,
    param,
    save_figures=bool(_config_value("bombcell", "save_plots", False)),
    return_figures=False,
)

missing = [path for path in _required_outputs() if not os.path.exists(path)]
if missing:
    raise FileNotFoundError(
        "Bombcell finished but required outputs are missing: "
        + ", ".join(missing)
    )

Path(str(snakemake.output.ready)).touch()
