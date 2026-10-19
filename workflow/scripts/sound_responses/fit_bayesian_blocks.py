"""Persist expensive raw segmentations independently of boundary selection."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(snakemake.scriptdir).resolve().parents[1]))
from utils.sound_responses.fits_session import fit_session

fit_session(snakemake.input.meta, snakemake.input.units, snakemake.output.h5,
            dict(snakemake.params.settings))
