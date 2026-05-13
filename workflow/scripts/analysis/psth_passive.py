import os
import sys
import json

import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D

# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.neurosuite import get_unit_names_sorted
from utils.psth import get_spike_counts
from utils.session import build_stimulus_catalog, detect_session_paradigm


def format_sound_label(sound_id, stimulus_catalog):
    entry = stimulus_catalog.get(str(sound_id), {})
    name = entry.get("name", f"sound_{sound_id}")
    freq = entry.get("freq", None)
    if freq is None:
        return name
    return f"{name} ({freq:g} Hz)"


def infer_stimulus_duration(sound_ids, stimulus_catalog):
    durations = []
    for sound_id in sound_ids:
        entry = stimulus_catalog.get(str(sound_id), {})
        duration = entry.get("duration", None)
        if duration is None:
            continue
        durations.append(float(duration))

    if not durations:
        return None
    return float(np.median(durations))


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


with h5py.File(snakemake.input[0], "r") as f:
    cfg = json.loads(f["processed"].attrs["parameters"])
    sound_events = np.array(f["processed"]["sound_events"])
    tl = np.array(f["processed"]["timeline"])

    paradigm = f["processed"].attrs.get("session_paradigm", detect_session_paradigm(cfg))
    if isinstance(paradigm, bytes):
        paradigm = paradigm.decode()

    stimulus_catalog_raw = f["processed"].attrs.get("stimulus_catalog", "")
    if isinstance(stimulus_catalog_raw, bytes):
        stimulus_catalog_raw = stimulus_catalog_raw.decode()
    stimulus_catalog = json.loads(stimulus_catalog_raw) if stimulus_catalog_raw else build_stimulus_catalog(cfg)

if paradigm != "passive":
    raise ValueError("psth_passive.py only supports passive sessions.")

spike_times = {}
with h5py.File(snakemake.input[1], "r") as f:
    unit_names = get_unit_names_sorted([name for name in f])
    for unit_name in unit_names:
        spike_times[unit_name] = np.sort(np.array(f[unit_name]["spike_times"]))

sound_ids = sound_events[:, 1].astype(np.int32)
event_ids = np.arange(len(sound_events), dtype=np.int32)
positive_sound_ids = sorted({int(sound_id) for sound_id in sound_ids if int(sound_id) > 0})

if not positive_sound_ids:
    for out_path in snakemake.output:
        write_placeholder_pdf(
            out_path,
            "Passive PSTH",
            "No positive sound IDs were found in processed/sound_events.",
        )
    raise SystemExit(0)

speed_thresh = cfg.get("position", {}).get("hd_update_speed", 0.04)
speed_ev = tl[sound_events[:, 2].astype(np.int32), 3]

event_groups = [
    ("Passive PSTH: all sounds", snakemake.output[0], event_ids),
    ("Passive PSTH: stationary sounds", snakemake.output[1], np.where(speed_ev < speed_thresh)[0]),
    ("Passive PSTH: running sounds", snakemake.output[2], np.where(speed_ev >= speed_thresh)[0]),
]

hw = snakemake.config["psth"]["micro"]["latency"]
bc = snakemake.config["psth"]["micro"]["bin_count"]
units_per_page = 12
cols = 3

color_map = plt.get_cmap("tab20")
sound_colors = {
    sound_id: color_map(i % color_map.N)
    for i, sound_id in enumerate(positive_sound_ids)
}

for title, out_path, idxs_group in event_groups:
    idxs_group = np.asarray(idxs_group, dtype=np.int32)
    event_idxs_by_sound = {}
    for sound_id in positive_sound_ids:
        idxs_for_sound = idxs_group[sound_ids[idxs_group] == sound_id]
        if len(idxs_for_sound) > 0:
            event_idxs_by_sound[sound_id] = idxs_for_sound

    plotted_sound_ids = list(event_idxs_by_sound.keys())
    if not plotted_sound_ids:
        write_placeholder_pdf(
            out_path,
            title,
            "No sound events were available for this movement state.",
        )
        continue

    legend_handles = [
        Line2D(
            [0],
            [0],
            color=sound_colors[sound_id],
            lw=2,
            label=f"{format_sound_label(sound_id, stimulus_catalog)} (n={len(event_idxs_by_sound[sound_id])})",
        )
        for sound_id in plotted_sound_ids
    ]
    stim_dur = infer_stimulus_duration(plotted_sound_ids, stimulus_catalog)

    with PdfPages(out_path) as pdf:
        for batch_start in range(0, len(unit_names), units_per_page):
            units_batch = unit_names[batch_start:batch_start + units_per_page]
            rows = int(np.ceil(len(units_batch) / cols))
            fig = plt.figure(figsize=(5 * cols, 3.6 * rows))

            for i, unit_name in enumerate(units_batch):
                ax = fig.add_subplot(rows, cols, i + 1)

                for sound_id in plotted_sound_ids:
                    event_times = sound_events[event_idxs_by_sound[sound_id], 0]
                    bins, psth = get_spike_counts(
                        spike_times[unit_name],
                        event_times,
                        hw=hw,
                        bin_count=bc,
                    )
                    bin_centers = bins[:-1] + (bins[1] - bins[0]) / 2.0
                    ax.plot(
                        bin_centers,
                        psth,
                        color=sound_colors[sound_id],
                        lw=1.6,
                        alpha=0.9,
                    )

                ax.axvline(0, color="black", ls="--", lw=1)
                if stim_dur is not None:
                    ax.axvspan(0, stim_dur, alpha=0.12, color="gray")
                ax.set_xlim(-hw, hw)
                ax.set_ylim(bottom=0)
                ax.set_title(unit_name, fontsize=11)
                ax.grid(alpha=0.2)
                if i % cols == 0:
                    ax.set_ylabel("Firing Rate, Hz", fontsize=10)
                if i >= len(units_batch) - cols:
                    ax.set_xlabel("Time, s", fontsize=10)

            fig.suptitle(title, fontsize=14)
            fig.legend(
                handles=legend_handles,
                loc="upper center",
                ncol=min(4, len(legend_handles)),
                fontsize=8,
                frameon=False,
            )
            fig.tight_layout(rect=(0, 0, 1, 0.93))
            pdf.savefig(fig)
            plt.close(fig)
