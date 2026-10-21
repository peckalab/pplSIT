"""Snakemake entry point for optional neural egocentric tuning."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(snakemake.scriptdir).resolve().parents[1]))
from utils.egocentric_tuning.session import run_session

run_session({name: str(getattr(snakemake.input, name)) for name in ('geometry', 'responses', 'units')},
            snakemake.output.h5, snakemake.output.csv, snakemake.output.comparisons,
            dict(snakemake.params.settings), str(snakemake.wildcards.session))
