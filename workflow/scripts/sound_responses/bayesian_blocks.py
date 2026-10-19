"""Snakemake entry point; numerical implementation lives in utils.sound_responses."""
from pathlib import Path
import sys

# Snakemake copies scripts to .snakemake/scripts; scriptdir retains the source location.
sys.path.insert(0, str(Path(snakemake.scriptdir).resolve().parents[1]))

from utils.sound_responses.session import select_session

select_session(snakemake.input.fits, snakemake.output.h5,
               snakemake.output.csv, dict(snakemake.params.settings))
