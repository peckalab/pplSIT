import json
import os
import sys

import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.aeps import AEP_metrics_lims
from utils.session import build_stimulus_catalog, detect_session_paradigm


def mean_ste(x):
    x = np.asarray(x, dtype=float)

    if x.ndim != 2 or x.shape[0] == 0:
        return np.full(x.shape[1] if x.ndim == 2 else 0, np.nan), np.full(x.shape[1] if x.ndim == 2 else 0, np.nan), 0

    mean = np.nanmean(x, axis=0)
    valid_trials = np.any(np.isfinite(x), axis=1)
    n_valid = int(np.sum(valid_trials))

    if n_valid == 0:
        ste = np.full(x.shape[1], np.nan)
    else:
        ste = np.nanstd(x[valid_trials], axis=0) / np.sqrt(n_valid)

    return mean, ste, n_valid


def _decode_attr(value):
    if isinstance(value, bytes):
        return value.decode()
    return value


def _load_meta(meta_path):
    with h5py.File(meta_path, "r") as f:
        tl = np.array(f["processed"]["timeline"])
        sound_events = np.array(f["processed"]["sound_events"])
        cfg = json.loads(f["processed"].attrs["parameters"])

        paradigm = _decode_attr(f["processed"].attrs.get("session_paradigm", detect_session_paradigm(cfg)))
        stimulus_catalog_raw = _decode_attr(f["processed"].attrs.get("stimulus_catalog", ""))

    if stimulus_catalog_raw:
        stimulus_catalog = json.loads(stimulus_catalog_raw)
    else:
        stimulus_catalog = build_stimulus_catalog(cfg)

    return tl, sound_events, cfg, paradigm, stimulus_catalog


def _format_sound_label(sound_id, stimulus_catalog):
    entry = stimulus_catalog.get(str(int(sound_id)), {})
    name = entry.get("name", f"sound_{int(sound_id)}")
    freq = entry.get("freq")
    if freq is None:
        return name
    return f"{name} ({freq:g} Hz)"


def _write_placeholder_pdf(pdf_path, title, message):
    fig = plt.figure(figsize=(8, 4))
    ax = fig.add_subplot(111)
    ax.axis("off")
    ax.text(0.5, 0.62, title, ha="center", va="center", fontsize=14)
    ax.text(0.5, 0.42, message, ha="center", va="center", fontsize=11)
    fig.tight_layout()
    with PdfPages(pdf_path) as pdf:
        pdf.savefig(fig)
    plt.close(fig)


def _plot_trace_group(ax, traces_by_condition, labels_by_condition, colors_by_condition, x_ms, area_name):
    plotted = 0
    for traces, label, color in zip(traces_by_condition, labels_by_condition, colors_by_condition):
        if traces.shape[0] == 0:
            continue

        mean, ste, n_valid = mean_ste(traces)
        if n_valid == 0 or mean.size == 0 or not np.any(np.isfinite(mean)):
            continue

        if np.any(np.isfinite(ste)):
            ax.fill_between(
                x_ms,
                0.2 * (mean - ste),
                0.2 * (mean + ste),
                color=color,
                alpha=0.3,
            )

        ax.plot(
            x_ms,
            0.2 * mean,
            color=color,
            lw=1.6,
            label=f"{label} (n={n_valid})",
        )
        plotted += 1

    ax.axhline(0, color="black", lw=1)
    ax.axvline(0, color="black", lw=1)
    ax.grid(alpha=0.25)

    if area_name in AEP_metrics_lims:
        for _, value in AEP_metrics_lims[area_name].items():
            ax.axvline(value[0], color="black", ls="--", lw=1)
            ax.axvline(value[1], color="black", ls="--", lw=1)

    if plotted == 0:
        ax.text(0.5, 0.5, "No valid trials", ha="center", va="center", transform=ax.transAxes)
    else:
        ax.legend(loc="upper right", fontsize=9)


tl, sound_events, cfg, paradigm, stimulus_catalog = _load_meta(snakemake.input[0])
if paradigm != "passive":
    raise ValueError("aeps_profiles_passive.py only supports passive sessions.")

aeps = {}
with h5py.File(snakemake.input[1], "r") as f:
    for area in sorted(f.keys()):
        aeps[area] = np.array(f[area]["avg_across_channels"])

if not aeps:
    _write_placeholder_pdf(snakemake.output[0], "Passive AEP profiles", "No area groups found in AEP file.")
    raise SystemExit(0)

sound_ids = sound_events[:, 1].astype(np.int32)
positive_sound_ids = sorted({int(sound_id) for sound_id in sound_ids if int(sound_id) > 0})
if not positive_sound_ids:
    _write_placeholder_pdf(
        snakemake.output[0],
        "Passive AEP profiles",
        "No positive sound IDs were found in processed/sound_events.",
    )
    raise SystemExit(0)

event_tl_idxs = sound_events[:, 2].astype(np.int32)
speed_ev = np.full(len(sound_events), np.nan, dtype=float)
valid_tl_mask = (event_tl_idxs >= 0) & (event_tl_idxs < len(tl))
speed_ev[valid_tl_mask] = tl[event_tl_idxs[valid_tl_mask], 3]

stationary_thresh = snakemake.config["lfp"]["baseline"]["stationary_thresh"]
stationary_mask = np.isfinite(speed_ev) & (speed_ev < stationary_thresh)
running_mask = np.isfinite(speed_ev) & (speed_ev >= stationary_thresh)

sound_pairs = [positive_sound_ids[i:i + 2] for i in range(0, len(positive_sound_ids), 2)]
aep_dur_ms = float(snakemake.config["lfp"]["aep_dur"]) * 1000.0

row_labels = [
    "All sounds in pair",
    "Stationary only",
    "Sound 1: stationary vs running",
    "Sound 2: stationary vs running",
]

pair_colors = ["tab:blue", "tab:orange"]
state_colors = ["navy", "orangered"]

with PdfPages(snakemake.output[0]) as pdf:
    for area, aeps_mx in aeps.items():
        cols = max(len(sound_pairs), 1)
        fig, axes = plt.subplots(4, cols, figsize=(5.2 * cols, 13), squeeze=False)

        for col, sound_pair in enumerate(sound_pairs):
            sound1 = sound_pair[0]
            sound2 = sound_pair[1] if len(sound_pair) > 1 else None

            label1 = _format_sound_label(sound1, stimulus_catalog)
            label2 = _format_sound_label(sound2, stimulus_catalog) if sound2 is not None else "No paired sound"

            pair_title = label1 if sound2 is None else f"{label1} / {label2}"
            axes[0, col].set_title(pair_title, fontsize=12)

            x_ms = np.linspace(0, aep_dur_ms, aeps_mx.shape[1])

            mask_sound1 = sound_ids == sound1
            mask_sound2 = sound_ids == sound2 if sound2 is not None else np.zeros(len(sound_ids), dtype=bool)

            _plot_trace_group(
                axes[0, col],
                [aeps_mx[mask_sound1], aeps_mx[mask_sound2]],
                [label1, label2],
                pair_colors,
                x_ms,
                area,
            )

            _plot_trace_group(
                axes[1, col],
                [aeps_mx[mask_sound1 & stationary_mask], aeps_mx[mask_sound2 & stationary_mask]],
                [f"{label1} sta", f"{label2} sta"],
                pair_colors,
                x_ms,
                area,
            )

            _plot_trace_group(
                axes[2, col],
                [aeps_mx[mask_sound1 & stationary_mask], aeps_mx[mask_sound1 & running_mask]],
                [f"{label1} sta", f"{label1} run"],
                state_colors,
                x_ms,
                area,
            )

            if sound2 is None:
                axes[3, col].axis("off")
                axes[3, col].text(0.5, 0.5, "No paired sound", ha="center", va="center", transform=axes[3, col].transAxes)
            else:
                _plot_trace_group(
                    axes[3, col],
                    [aeps_mx[mask_sound2 & stationary_mask], aeps_mx[mask_sound2 & running_mask]],
                    [f"{label2} sta", f"{label2} run"],
                    state_colors,
                    x_ms,
                    area,
                )

            for row in range(4):
                ax = axes[row, col]
                if row < 3 or sound2 is not None:
                    ax.set_xlabel("Time from stim. onset, ms", fontsize=10)

            if col == 0:
                for row, row_label in enumerate(row_labels):
                    axes[row, col].set_ylabel(f"{row_label}\nLFP, $\\mu$V", fontsize=10)

        fig.suptitle(f"{area} passive AEP profiles", fontsize=16)
        fig.tight_layout(rect=[0, 0, 1, 0.97])
        pdf.savefig(fig)
        plt.close(fig)
