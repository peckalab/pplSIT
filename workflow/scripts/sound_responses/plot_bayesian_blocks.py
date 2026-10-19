"""Plot saved results without rerunning Bayesian-block estimation."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(snakemake.scriptdir).resolve().parents[1]))

from utils.sound_responses.plots import plot_results

plot_results(snakemake.input.h5, snakemake.output.pdf)
