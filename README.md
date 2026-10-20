# pplSIT
Pipeline for SIT data postprocessing


Example execution (from the 'workflow' folder):
snakemake --configfile="../config/andrey_pd.yaml" --cores=1

MoSeq heading coordinates and regeneration of earlier synchronized outputs
are described in [the synchronization notes](workflow/scripts/README_moseq_sync.md).

The optional [egocentric geometry stage](workflow/scripts/egocentric_geometry/README.md)
prepares head/body-relative speaker bearings and event-coverage plots for spatial tuning.
