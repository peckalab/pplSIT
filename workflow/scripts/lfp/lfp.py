import os
import h5py
import numpy as np
from scipy import signal


def get_staged_stream_dat(snakemake):
    """
    Resolve the staged .dat file for the current session/stream.

    Expected staged location:
        <src_path>/<animal>/<session>/ephys/<stream>/*.dat
    """
    animal = snakemake.wildcards.animal
    session = snakemake.wildcards.session
    stream = snakemake.wildcards.stream

    ephys_stream_dir = os.path.join(
        snakemake.config["src_path"], animal, session, "ephys", stream
    )

    if not os.path.isdir(ephys_stream_dir):
        raise FileNotFoundError(
            f"Staged stream directory does not exist: {ephys_stream_dir}"
        )

    dat_files = sorted(
        os.path.join(ephys_stream_dir, f)
        for f in os.listdir(ephys_stream_dir)
        if f.endswith(".dat")
    )

    if len(dat_files) == 0:
        raise FileNotFoundError(
            f"No .dat files found in staged stream directory: {ephys_stream_dir}"
        )

    if len(dat_files) > 1:
        raise ValueError(
            f"Expected exactly one .dat file in {ephys_stream_dir}, found {len(dat_files)}:\n"
            + "\n".join(dat_files)
        )

    return dat_files[0]


# -----------------------
# config
# -----------------------
dtype = np.int16  # TODO make configurable?
fs = snakemake.config["lfp"]["source_rate"]
lfp_fs = snakemake.config["lfp"]["target_rate"]
n_channels = snakemake.config["lfp"]["n_channels"]
chunk_size_sec = snakemake.config["lfp"]["chunk_size_sec"]

channels_to_extract = None
if channels_to_extract is None:
    channels_to_extract = list(range(n_channels))

ds_factor = fs // lfp_fs
assert fs % lfp_fs == 0, "fs must be divisible by lfp_fs"

# Design lowpass filter (~0.8 * Nyquist of new fs)
sos = signal.butter(4, 0.8 * (lfp_fs / 2) / (fs / 2), btype="low", output="sos")

# -----------------------
# resolve input .dat
# -----------------------
dat_path = get_staged_stream_dat(snakemake)
print(f"LFP source .dat: {dat_path}")

# Determine total number of samples
n_samples_total = int(np.memmap(dat_path, dtype=dtype).size / n_channels)
n_samples_lfp = n_samples_total // ds_factor
n_channels_out = len(channels_to_extract)
print_perc = n_samples_total / 10

# -----------------------
# create output HDF5
# -----------------------
os.makedirs(os.path.dirname(snakemake.output[0]), exist_ok=True)

with h5py.File(snakemake.output[0], "w") as h5f:
    dset = h5f.create_dataset(
        "lfp",
        shape=(n_channels_out, n_samples_lfp),
        dtype="float32"
    )

    samples_per_chunk = int(chunk_size_sec * fs)
    lfp_idx = 0

    with open(dat_path, "rb") as f:
        for start in range(0, n_samples_total, samples_per_chunk):
            stop = min(start + samples_per_chunk, n_samples_total)
            n_samples = stop - start

            # Read chunk
            f.seek(start * n_channels * np.dtype(dtype).itemsize)
            data = np.fromfile(f, dtype=dtype, count=n_samples * n_channels)
            data = data.reshape((n_samples, n_channels)).T
            data = data[channels_to_extract, :]

            # Filter and downsample
            data_filt = signal.sosfiltfilt(sos, data, axis=1)
            data_ds = data_filt[:, ::ds_factor]

            # Store in HDF5
            available = dset.shape[1] - lfp_idx
            write_len = min(available, data_ds.shape[1])
            dset[:, lfp_idx:lfp_idx + write_len] = data_ds[:, :write_len]
            lfp_idx += write_len

            if start > print_perc:
                print(f"LFP: {start / n_samples_total * 100:.1f}% done")
                print_perc += n_samples_total / 10