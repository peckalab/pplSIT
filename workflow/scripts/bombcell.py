import json
import os
import fcntl
import glob
import shutil
import sys
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

sys.path = [
    path for path in sys.path
    if not (Path(path or ".") / "bombcell.py").exists()
]
import bombcell


def _config_value(section, key, default=None):
    return snakemake.config.get(section, {}).get(key, default)


def _acquire_io_heavy_lock():
    lock_config = snakemake.config.get("io_heavy_lock", {})
    if not lock_config.get("enabled", False):
        return None

    lock_path = str(lock_config.get("path", "/tmp/pplSIT_io_heavy.lock"))
    os.makedirs(os.path.dirname(lock_path), exist_ok=True)
    lock_handle = open(lock_path, "w")
    print(f"Waiting for shared heavy-I/O lock: {lock_path}", flush=True)
    fcntl.flock(lock_handle, fcntl.LOCK_EX)
    print(f"Acquired shared heavy-I/O lock: {lock_path}", flush=True)
    return lock_handle


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


def _single_match(pattern, description):
    matches = sorted(glob.glob(pattern))
    if len(matches) != 1:
        raise FileNotFoundError(
            f"Expected exactly one {description} matching {pattern}, found: {matches}"
        )
    return matches[0]


def _require_files(paths, description):
    missing = [path for path in paths if not os.path.exists(path)]
    if missing:
        raise FileNotFoundError(
            f"Missing {description}: " + ", ".join(missing)
        )


def _raw_waveform_files(save_path):
    return sorted((Path(save_path) / "RawWaveforms").glob("Unit*_RawSpikes.npy"))


def _write_raw_waveform_marker(raw_marker, raw_waveform_files, param):
    raw_marker = Path(raw_marker)
    raw_marker.parent.mkdir(parents=True, exist_ok=True)
    with raw_marker.open("w") as f:
        json.dump({
            "unit_match_waveforms": bool(_config_value("bombcell", "unit_match_waveforms", False)),
            "n_raw_waveform_files": len(raw_waveform_files),
            "nRawSpikesToExtract": int(param["nRawSpikesToExtract"]),
            "saveMultipleRaw": bool(param.get("saveMultipleRaw", False)),
            "spike_width": int(param["spike_width"]),
        }, f, indent=2)
        f.write("\n")


def _raw_waveforms_match_config(save_path):
    raw_marker = Path(str(snakemake.output.raw_waveforms_ready))
    if not raw_marker.exists():
        return False

    try:
        with raw_marker.open("r") as f:
            marker = json.load(f)
    except (OSError, json.JSONDecodeError):
        return False

    expected_spike_width = int(_config_value("unit_match", "spike_width", marker.get("spike_width", 0)))
    if int(marker.get("spike_width", -1)) != expected_spike_width:
        return False

    expected_sample_amount = int(
        _config_value("unit_match", "sample_amount", marker.get("nRawSpikesToExtract", 0))
    )
    if int(marker.get("nRawSpikesToExtract", -1)) != expected_sample_amount:
        return False

    if not bool(marker.get("saveMultipleRaw", False)):
        return False

    raw_waveform_files = _raw_waveform_files(save_path)
    if not raw_waveform_files:
        return False

    marker_n_files = marker.get("n_raw_waveform_files")
    if marker_n_files is not None and int(marker_n_files) != len(raw_waveform_files):
        return False

    return True


def _outputs_complete(save_path):
    required = _required_outputs() + [str(snakemake.output.ready)]
    if not all(os.path.exists(path) for path in required):
        return False

    if bool(_config_value("bombcell", "unit_match_waveforms", False)):
        if bool(_config_value("bombcell", "force_unit_match_waveform_reextract", False)):
            return False
        if not _raw_waveforms_match_config(save_path):
            return False

    return True


def _stage_raw_file(raw_file):
    stage_config = _config_value("bombcell", "stage_raw_dat", {}) or {}
    if not bool(stage_config.get("enabled", False)):
        return raw_file, None

    stage_root = Path(str(stage_config.get("root", "/local/scratch/bengala/pplSIT/bombcell_stage")))
    cleanup = bool(stage_config.get("cleanup", True))
    min_free_gb = float(stage_config.get("min_free_gb", 20))

    stage_name = "__".join([
        str(snakemake.wildcards.animal),
        str(snakemake.wildcards.session),
        str(snakemake.wildcards.stream),
    ])
    stage_dir = stage_root / stage_name
    stage_dir.mkdir(parents=True, exist_ok=True)

    raw_file = Path(raw_file)
    staged_raw = stage_dir / raw_file.name
    raw_size = raw_file.stat().st_size
    free_bytes = shutil.disk_usage(stage_dir).free
    required_free = raw_size + int(min_free_gb * 1024 ** 3)
    if not (staged_raw.exists() and staged_raw.stat().st_size == raw_size) and free_bytes < required_free:
        raise OSError(
            f"Not enough free space to stage {raw_file} into {stage_dir}. "
            f"Need at least {(required_free / 1024 ** 3):.1f} GiB, "
            f"found {(free_bytes / 1024 ** 3):.1f} GiB."
        )

    if staged_raw.exists() and staged_raw.stat().st_size == raw_size:
        print(f"Reusing staged raw file: {staged_raw}", flush=True)
    else:
        tmp_raw = staged_raw.with_suffix(staged_raw.suffix + ".part")
        if tmp_raw.exists():
            tmp_raw.unlink()
        print(f"Staging raw file to local scratch: {raw_file} -> {staged_raw}", flush=True)
        shutil.copy2(raw_file, tmp_raw)
        os.replace(tmp_raw, staged_raw)

    return str(staged_raw), stage_dir if cleanup else None


ks_dir = os.path.dirname(str(snakemake.input.st))
save_path = os.path.dirname(str(snakemake.output.ready))
settings_path = os.path.join(ks_dir, "settings.json")
io_heavy_lock_handle = _acquire_io_heavy_lock()

os.makedirs(save_path, exist_ok=True)

if _outputs_complete(save_path):
    print(
        "BombCell outputs and UnitMatch raw waveforms already exist; "
        "skipping extraction/QC for this stream.",
        flush=True,
    )
    raise SystemExit(0)

raw_file = _single_match(os.path.join(ks_dir, "*.dat"), "Kilosort stream .dat file")

_require_files(
    [
        str(snakemake.input.st),
        str(snakemake.input.sc),
        str(snakemake.input.templates),
        os.path.join(ks_dir, "spike_templates.npy"),
        os.path.join(ks_dir, "amplitudes.npy"),
        os.path.join(ks_dir, "whitening_mat_inv.npy"),
        os.path.join(ks_dir, "channel_positions.npy"),
        settings_path,
        raw_file,
    ],
    "Bombcell Kilosort inputs",
)

raw_file, stage_dir_to_cleanup = _stage_raw_file(raw_file)

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

if bool(_config_value("bombcell", "unit_match_waveforms", False)):
    raw_waveforms_dir = Path(save_path) / "RawWaveforms"
    raw_waveforms_complete = _raw_waveforms_match_config(save_path)
    force_reextract = bool(_config_value("bombcell", "force_unit_match_waveform_reextract", False))
    if (force_reextract or not raw_waveforms_complete) and raw_waveforms_dir.exists():
        shutil.rmtree(raw_waveforms_dir)
    param["extractRaw"] = True
    param["reextractRaw"] = force_reextract or not raw_waveforms_complete
    param["nRawSpikesToExtract"] = int(
        _config_value("unit_match", "sample_amount", param["nRawSpikesToExtract"])
    )
    param["saveMultipleRaw"] = True
    param["detrendForUnitMatch"] = False
    param["spike_width"] = int(_config_value("unit_match", "spike_width", param["spike_width"]))
    if param["spike_width"] == 82:
        param["waveformBaselineNoiseWindow"] = 20
    param["decompress_data"] = bool(_config_value("bombcell", "decompress_data", param["decompress_data"]))

try:
    bombcell.run_bombcell(
        ks_dir,
        save_path,
        param,
        save_figures=bool(_config_value("bombcell", "save_plots", False)),
        return_figures=False,
    )
finally:
    if stage_dir_to_cleanup is not None:
        print(f"Cleaning staged raw file directory: {stage_dir_to_cleanup}", flush=True)
        shutil.rmtree(stage_dir_to_cleanup, ignore_errors=True)

missing = [path for path in _required_outputs() if not os.path.exists(path)]
if missing:
    raise FileNotFoundError(
        "Bombcell finished but required outputs are missing: "
        + ", ".join(missing)
    )

Path(str(snakemake.output.ready)).touch()

raw_waveform_files = _raw_waveform_files(save_path)
_write_raw_waveform_marker(str(snakemake.output.raw_waveforms_ready), raw_waveform_files, param)

if bool(_config_value("bombcell", "unit_match_waveforms", False)) and not raw_waveform_files:
    raise FileNotFoundError(
        "BombCell unit_match_waveforms=True but no RawWaveforms/Unit*_RawSpikes.npy files were created."
    )
