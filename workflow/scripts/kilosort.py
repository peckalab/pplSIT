import os, json

from kilosort import run_kilosort, DEFAULT_SETTINGS
from kilosort.io import load_probe


settings = DEFAULT_SETTINGS

# load settings
with open(snakemake.input[0]) as json_file:  
    settings_local = json.load(json_file)

settings.update(settings_local)
settings['data_dir'] = os.path.dirname(snakemake.input[2])

results_dir = os.path.dirname(snakemake.input[2])

# load probe configuration
probe = load_probe(snakemake.input[1])

save_p = snakemake.config['kilosort']['save_preprocessed']

# run kilosort
ops, st, clu, tF, Wall, similar_templates, is_ref, est_contam_rate, kept_spikes = \
    run_kilosort(settings=settings, probe=probe, results_dir=results_dir, save_preprocessed_copy=save_p)

# save configuration
with open(os.path.join(results_dir, 'settings.json'), 'w') as f:
    f.write(json.dumps(settings, indent=2))