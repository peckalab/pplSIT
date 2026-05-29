import json
import os
import sys

import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.neurosuite import get_unit_names_sorted
from utils.psth import get_spike_counts
from utils.session import build_stimulus_catalog, detect_session_paradigm


def cfg_get(name, default):
    return snakemake.config.get("strf_passive", {}).get(name, default)


def decode_attr(value):
    if isinstance(value, bytes):
        return value.decode()
    return value


def write_placeholder_pdf(pdf_path, title, message):
    fig = plt.figure(figsize=(8, 4))
    ax = fig.add_subplot(111)
    ax.axis("off")
    ax.text(0.5, 0.62, title, ha="center", va="center", fontsize=14)
    ax.text(0.5, 0.42, message, ha="center", va="center", fontsize=11)
    fig.tight_layout()
    with PdfPages(pdf_path) as pdf:
        pdf.savefig(fig)
    plt.close(fig)


def passive_stimuli(sound_events, stimulus_catalog):
    entries = []
    sound_ids = sorted({int(x) for x in sound_events[:, 1] if int(x) > 0})

    for sound_id in sound_ids:
        entry = stimulus_catalog.get(str(sound_id), {})
        if "freq" not in entry:
            continue
        entries.append(
            {
                "sound_id": sound_id,
                "freq": float(entry["freq"]),
                "name": entry.get("name", f"sound_{sound_id}"),
            }
        )

    return sorted(entries, key=lambda item: item["freq"])


def frequency_edges(freqs):
    freqs = np.asarray(freqs, dtype=np.float64)
    if len(freqs) == 1:
        return np.array([freqs[0] / np.sqrt(2.0), freqs[0] * np.sqrt(2.0)])

    log_freqs = np.log(freqs)
    mids = (log_freqs[:-1] + log_freqs[1:]) / 2.0
    first = log_freqs[0] - (mids[0] - log_freqs[0])
    last = log_freqs[-1] + (log_freqs[-1] - mids[-1])
    return np.exp(np.concatenate([[first], mids, [last]]))


def baseline_subtract(matrix, bins):
    baseline_mask = bins[:-1] < 0
    if not np.any(baseline_mask):
        return matrix.copy()
    baseline = np.nanmean(matrix[:, baseline_mask], axis=1, keepdims=True)
    return matrix - baseline


def robust_vlim(matrix):
    values = matrix[np.isfinite(matrix)]
    if values.size == 0:
        return 1.0
    vmax = np.nanpercentile(np.abs(values), 98)
    return float(vmax) if vmax > 0 else 1.0


with h5py.File(snakemake.input.meta, "r") as f:
    cfg = json.loads(f["processed"].attrs["parameters"])
    sound_events = np.array(f["processed"]["sound_events"])
    timeline = np.array(f["processed"]["timeline"])

    paradigm = decode_attr(f["processed"].attrs.get("session_paradigm", detect_session_paradigm(cfg)))
    stimulus_catalog_raw = decode_attr(f["processed"].attrs.get("stimulus_catalog", ""))
    stimulus_catalog = json.loads(stimulus_catalog_raw) if stimulus_catalog_raw else build_stimulus_catalog(cfg)

if paradigm != "passive":
    raise ValueError("strf_passive.py only supports passive sessions.")

with h5py.File(snakemake.input.units, "r") as f:
    unit_names = get_unit_names_sorted([name for name in f])
    spike_times = {
        unit_name: np.sort(np.array(f[unit_name]["spike_times"]))
        for unit_name in unit_names
    }

stimuli = passive_stimuli(sound_events, stimulus_catalog)
if not stimuli:
    os.makedirs(os.path.dirname(snakemake.output.h5), exist_ok=True)
    with h5py.File(snakemake.output.h5, "w") as f:
        f.attrs["message"] = "No positive passive sound IDs with frequencies were found."
    write_placeholder_pdf(
        snakemake.output.pdf,
        "Passive STRF-like maps",
        "No positive passive sound IDs with frequencies were found.",
    )
    raise SystemExit(0)

pre = float(cfg_get("pre", 0.05))
post = float(cfg_get("post", 0.25))
bin_count = int(cfg_get("bin_count", 61))
min_events = int(cfg_get("min_events", 3))
units_per_page = int(cfg_get("units_per_page", 8))
speed_thresh = float(cfg_get("speed_thresh", cfg.get("position", {}).get("hd_update_speed", 0.04)))

hw = max(pre, post)
bins_full = np.linspace(-hw, hw, bin_count)
bin_mask = (bins_full[:-1] >= -pre) & (bins_full[1:] <= post)
bins = np.concatenate([bins_full[:-1][bin_mask], [bins_full[1:][bin_mask][-1]]])
bin_centers = bins[:-1] + np.diff(bins) / 2.0

sound_ids = sound_events[:, 1].astype(np.int32)
speed_at_event = timeline[sound_events[:, 2].astype(np.int32), 3]
all_event_idxs = np.arange(len(sound_events), dtype=np.int32)

event_groups = [
    ("all", all_event_idxs),
    ("stationary", np.where(speed_at_event < speed_thresh)[0]),
    ("running", np.where(speed_at_event >= speed_thresh)[0]),
]

freqs = np.array([item["freq"] for item in stimuli], dtype=np.float64)
stimulus_ids = np.array([item["sound_id"] for item in stimuli], dtype=np.int32)
stimulus_names = [item["name"] for item in stimuli]

os.makedirs(os.path.dirname(snakemake.output.h5), exist_ok=True)
with h5py.File(snakemake.output.h5, "w") as f:
    f.attrs["description"] = (
        "Passive STRF-like frequency-by-latency maps from pure-tone PSTHs. "
        "These are not TORC/dynamic-ripple STRFs."
    )
    f.attrs["pre"] = pre
    f.attrs["post"] = post
    f.attrs["bin_count"] = bin_count
    f.attrs["min_events"] = min_events
    f.attrs["speed_thresh"] = speed_thresh
    f.create_dataset("bins", data=bins)
    f.create_dataset("bin_centers", data=bin_centers)
    f.create_dataset("frequencies", data=freqs)
    f.create_dataset("sound_ids", data=stimulus_ids)
    f.attrs["sound_names"] = json.dumps(stimulus_names)

    for group_name, idxs_group in event_groups:
        group = f.create_group(group_name)
        event_counts = np.zeros(len(stimuli), dtype=np.int32)

        for unit_name in unit_names:
            raw_matrix = np.full((len(stimuli), len(bin_centers)), np.nan, dtype=np.float64)

            for row_idx, stimulus in enumerate(stimuli):
                event_idxs = idxs_group[sound_ids[idxs_group] == stimulus["sound_id"]]
                event_counts[row_idx] = len(event_idxs)
                if len(event_idxs) < min_events:
                    continue

                event_times = sound_events[event_idxs, 0]
                _, psth_full = get_spike_counts(
                    spike_times[unit_name],
                    event_times,
                    hw=hw,
                    bin_count=bin_count,
                )
                raw_matrix[row_idx] = psth_full[bin_mask]

            unit_group = group.create_group(unit_name)
            unit_group.create_dataset("rate_hz", data=raw_matrix)
            unit_group.create_dataset("baseline_subtracted_hz", data=baseline_subtract(raw_matrix, bins))

        group.create_dataset("event_counts", data=event_counts)

freq_edges = frequency_edges(freqs)

with h5py.File(snakemake.output.h5, "r") as f_in, PdfPages(snakemake.output.pdf) as pdf:
    for group_name, _ in event_groups:
        if group_name not in f_in:
            continue

        group = f_in[group_name]
        event_counts = np.array(group["event_counts"])
        valid_freqs = event_counts >= min_events

        for batch_start in range(0, len(unit_names), units_per_page):
            units_batch = unit_names[batch_start:batch_start + units_per_page]
            cols = 2
            rows = int(np.ceil(len(units_batch) / cols))
            fig = plt.figure(figsize=(6.0 * cols, 3.8 * rows))

            for i, unit_name in enumerate(units_batch):
                ax = fig.add_subplot(rows, cols, i + 1)
                matrix = np.array(group[unit_name]["baseline_subtracted_hz"])
                matrix_plot = np.ma.masked_invalid(matrix)
                vlim = robust_vlim(matrix)

                pcm = ax.pcolormesh(
                    bins,
                    freq_edges,
                    matrix_plot,
                    cmap="RdBu_r",
                    vmin=-vlim,
                    vmax=vlim,
                    shading="auto",
                )
                ax.axvline(0, color="black", lw=1, ls="--")
                ax.set_yscale("log")
                ax.set_yticks(freqs)
                ax.set_yticklabels([f"{freq:g}" for freq in freqs], fontsize=8)
                ax.set_xlim(-pre, post)
                ax.set_title(unit_name, fontsize=11)
                ax.set_xlabel("Latency, s", fontsize=9)
                ax.set_ylabel("Frequency, Hz", fontsize=9)

                if not np.any(valid_freqs):
                    ax.text(
                        0.5,
                        0.5,
                        "No frequencies reached min_events",
                        ha="center",
                        va="center",
                        transform=ax.transAxes,
                        fontsize=9,
                    )

                cbar = fig.colorbar(pcm, ax=ax, fraction=0.046, pad=0.04)
                cbar.set_label("Delta firing rate, Hz", fontsize=8)
                cbar.ax.tick_params(labelsize=8)

            fig.suptitle(
                f"Passive STRF-like maps: {group_name} events",
                fontsize=14,
            )
            fig.tight_layout(rect=(0, 0, 1, 0.95))
            pdf.savefig(fig)
            plt.close(fig)
