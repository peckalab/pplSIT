"""Snakemake entry point for optional egocentric geometry."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(snakemake.scriptdir).resolve().parents[1]))
from utils.egocentric_geometry.session import run_session

inputs = {name: str(getattr(snakemake.input, name))
          for name in ('meta', 'dlc')}
inputs['moseq'] = str(snakemake.input.moseq[0]) if snakemake.input.moseq else None
run_session(inputs, snakemake.output.h5, snakemake.output.csv,
            dict(snakemake.params.settings), str(snakemake.wildcards.session))
