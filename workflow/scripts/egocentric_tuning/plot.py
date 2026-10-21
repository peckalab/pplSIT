"""Snakemake entry point for neural tuning plots."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(snakemake.scriptdir).resolve().parents[1]))
from utils.egocentric_tuning.plots import plot_tuning

plot_tuning(snakemake.input.h5, snakemake.output.pdf, snakemake.output.preview)
