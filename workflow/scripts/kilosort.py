import os
import sys
import json
import glob
import torch
import numpy as np

from kilosort import run_kilosort, DEFAULT_SETTINGS
from kilosort.io import load_probe

# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.kilosort import infer_n_chan_bin


def get_stream_kilosort_dir(snakemake):
    """
    Return processed/<animal>/<session>/kilosort/<stream>
    """
    return os.path.join(
        snakemake.config["dst_path"],
        snakemake.wildcards.animal,
        snakemake.wildcards.session,
        "kilosort",
        snakemake.wildcards.stream,
    )


def resolve_stream_inputs(snakemake):
    """
    Resolve settings.json, probe.json, and the single .dat file
    from the current stream folder.
    """
    stream_dir = get_stream_kilosort_dir(snakemake)

    if not os.path.isdir(stream_dir):
        raise FileNotFoundError(f"Kilosort stream directory not found: {stream_dir}")

    settings_path = os.path.join(stream_dir, "settings.json")
    probe_path = os.path.join(stream_dir, "probe.json")

    if not os.path.exists(settings_path):
        raise FileNotFoundError(f"settings.json not found: {settings_path}")

    if not os.path.exists(probe_path):
        raise FileNotFoundError(f"probe.json not found: {probe_path}")

    dat_files = sorted(glob.glob(os.path.join(stream_dir, "*.dat")))
    if len(dat_files) == 0:
        raise FileNotFoundError(f"No .dat found in {stream_dir}")
    if len(dat_files) > 1:
        raise ValueError(f"Expected exactly one .dat in {stream_dir}, found: {dat_files}")

    dat_path = dat_files[0]
    return settings_path, probe_path, dat_path, stream_dir


settings = DEFAULT_SETTINGS.copy()

# --- resolve paths from stream folder ---
settings_path, probe_path, dat_path, stream_dir = resolve_stream_inputs(snakemake)

# --- outputs ---
# Put results next to spike_times.npy (stream folder)
results_dir = os.path.dirname(snakemake.output["st"])

# load settings
with open(settings_path) as f:
    settings_local = json.load(f)

settings.update(settings_local)

# Kilosort expects data_dir pointing to folder containing .dat
settings["data_dir"] = os.path.dirname(dat_path)

n_chan_bin = infer_n_chan_bin(dat_path, base_n_chan=384, dtype_bytes=2)
settings["n_chan_bin"] = n_chan_bin

# load probe configuration
probe = load_probe(probe_path)

save_p = snakemake.config["kilosort"].get("save_preprocessed", False)

# select CUDA device with most free memory (fallback to configured device)
def pick_best_cuda_device(fallback_device=0, min_free_gb=8.0, max_utilization=20.0):
    if not torch.cuda.is_available():
        print("CUDA not available")
        return fallback_device

    try:
        import subprocess

        print("CUDA_VISIBLE_DEVICES =", os.environ.get("CUDA_VISIBLE_DEVICES"))

        cmd = [
            "nvidia-smi",
            "--query-gpu=index,memory.free,memory.total,utilization.gpu",
            "--format=csv,noheader,nounits",
        ]
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)

        gpu_rows = []
        for line in result.stdout.strip().splitlines():
            idx_s, free_s, total_s, util_s = [x.strip() for x in line.split(",")]
            g = {
                "idx": int(idx_s),
                "free_mib": float(free_s),
                "total_mib": float(total_s),
                "util": float(util_s),
            }
            g["free_frac"] = g["free_mib"] / max(g["total_mib"], 1.0)
            gpu_rows.append(g)

        print("GPU rows from nvidia-smi:")
        for g in gpu_rows:
            print(g)

        min_free_mib = min_free_gb * 1024.0
        preferred = [
            g for g in gpu_rows
            if g["util"] <= max_utilization and g["free_mib"] >= min_free_mib
        ]

        print("Preferred GPUs:")
        for g in preferred:
            print(g)

        if preferred:
            best = max(preferred, key=lambda g: (g["free_mib"], -g["util"]))
            print("Selected from preferred:", best)
            return best["idx"]

        for g in gpu_rows:
            g["score"] = (100.0 - g["util"]) * 2 + 50.0 * g["free_frac"]

        print("Scored GPUs:")
        for g in gpu_rows:
            print(g)

        best = max(gpu_rows, key=lambda g: g["score"])
        print("Selected from fallback:", best)
        return best["idx"]

    except Exception as e:
        print(f"GPU auto-selection failed: {e}")
        return fallback_device

# select CUDA device using utilization + memory (fallback to configured device)
best_dev_id = snakemake.config["kilosort"].get("cuda_device", 0)
best_dev_id = pick_best_cuda_device(
    fallback_device=best_dev_id,
    min_free_gb=8.0,
    max_utilization=20.0,
)

print(f"USING CUDA DEVICE: {best_dev_id}")
print(f"KILOSORT STREAM DIR: {stream_dir}")
print(f"KILOSORT DAT: {dat_path}")
print(f"KILOSORT SETTINGS: {settings_path}")
print(f"KILOSORT PROBE: {probe_path}")

assert probe["n_chan"] == 384 or getattr(probe, "n_chan", None) == 384
if n_chan_bin == 385:
    print("Detected 385 channels in binary; assuming last channel is sync and excluded by chanMap.")

# run kilosort
ops, st, clu, tF, Wall, similar_templates, is_ref, est_contam_rate, kept_spikes = run_kilosort(
    settings=settings,
    probe=probe,
    results_dir=results_dir,
    save_preprocessed_copy=save_p,
    device=torch.device(best_dev_id) if torch.cuda.is_available() else torch.device("cpu"),
)

# save configuration actually used
with open(os.path.join(results_dir, "settings_used.json"), "w") as f:
    json.dump(settings, f, indent=2)


"""
* Resulting files in the dataset:
    * _First a note about many of these files:_ Many of the files are `.npy` format which is a simple file format 
        that can be natively read in python and read in matlab via the [npy-matlab repository](https://github.com/kwikteam/npy-matlab).
    * A raw data file, with any filename. This file should be "flat binary" format, meaning that the data values 
        corresponding to the voltage traces can are just the literal bytes in the file with no additional formatting 
        (header data is allowed at the beginning of the file, but will not be used)
    * `params.py` - text file that specifies:
        * `dat_path` - location of raw data file
        * `n_channels_dat` - _total_ number of rows in the data file (not just those that have your neural data on them. 
            This is for loading the file)
        * `dtype` - data type to read, e.g. 'int16'
        * `offset` - number of bytes at the beginning of the file to skip
        * `sample_rate` - in Hz
        * `hp_filtered` - True/False, whether the data have already been filtered
    * `amplitudes.npy` - `[nSpikes, ] double` vector with the amplitude scaling factor that was applied to the template when extracting that spike
    * `channel_map.npy` - `[nChannels, ] int32` vector with the channel map, i.e. which row of the data file to look in for the channel in question
    * `channel_positions.npy` - `[nChannels, 2] double` matrix with each row giving the x and y coordinates of that channel. 
        Together with the channel map, this determines how waveforms will be plotted in WaveformView (see below).
    * `pc_features.npy` - `[nSpikes, nFeaturesPerChannel, nPCFeatures] single` matrix giving the PC values for each spike. 
        The channels that those features came from are specified in pc_features_ind.npy. E.g. the value at 
        `pc_features[123, 1, 5]` is the projection of the 123rd spike onto the 1st PC on the channel given by `pc_feature_ind[5]`.
    * `pc_feature_ind.npy` - `[nTemplates, nPCFeatures] uint32` matrix specifying which pcFeatures are included in the pc_features matrix.
    * `similar_templates.npy` - `[nTemplates, nTemplates] single` matrix giving the similarity score (larger is more similar) 
        between each pair of templates
    * `spike_templates.npy` - `[nSpikes, ] uint32` vector specifying the identity of the template that was used to extract each spike
    * `spike_times.npy` - `[nSpikes, ] uint64` vector giving the spike time of each spike in **samples**. To convert to seconds, 
        divide by sample_rate from params.py.
    * `template_features.npy` - `[nSpikes, nTempFeatures] single` matrix giving the magnitude of the projection 
        of each spike onto nTempFeatures other features. Which other features is specified in `template_feature_ind.npy`
    * `template_feature_ind.npy` - `[nTemplates, nTempFeatures] uint32` matrix specifying which templateFeatures are 
        included in the template_features matrix.
    * `templates.npy` - `[nTemplates, nTimePoints, nTempChannels] single` matrix giving the template shapes on the 
        channels given in `templates_ind.npy`
    * `templates_ind.npy` - `[nTemplates, nTempChannels] double` matrix specifying the channels on which each template is defined. 
        In the case of Kilosort templates_ind is just the integers from 0 to nChannels-1, since templates are defined on all channels.
    * `whitening_mat.npy` - `[nChannels, nChannels] double` whitening matrix applied to the data during automatic spike sorting
    * `whitening_mat_inv.npy` - `[nChannels, nChannels] double`, the inverse of the whitening matrix.
    * `spike_clusters.npy` - `[nSpikes, ] int32` vector giving the cluster identity of each spike. This file is optional 
        and if not provided will be automatically created the first time you run the template gui, taking the same values as 
        spike_templates.npy until you do any merging or splitting.
    * `cluster_groups.csv` - comma-separated value text file giving the "cluster group" of each cluster (0=noise, 1=MUA, 2=Good, 3=unsorted)
"""