"""Render the familiar BGR/TGT bar plot with saved response timing."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(snakemake.scriptdir).resolve().parents[1]))
from utils.sound_responses.psth_overlay import plot_psth_micro_overlay

optional = list(snakemake.input.phase_lock)
plot_psth_micro_overlay(snakemake.input.meta, snakemake.input.units,
                        snakemake.input.responses, snakemake.output.pdf,
                        dict(snakemake.params.micro),
                        phase_lock_path=optional[0] if optional else None,
                        preview_png=snakemake.output.preview)
