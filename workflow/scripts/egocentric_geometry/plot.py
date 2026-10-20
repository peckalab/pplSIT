"""Snakemake entry point for geometry coverage plots."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(snakemake.scriptdir).resolve().parents[1]))
from utils.egocentric_geometry.plots import plot_geometry

plot_geometry(snakemake.input.h5, snakemake.output.pdf, snakemake.output.preview)
