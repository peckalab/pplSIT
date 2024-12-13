import os, sys, json
import numpy as np


# main function to create all inputs and wildcards
def get_clu_input_names(w):
    # TODO: this should be done in a more elegant way
    if os.path.isfile(os.path.join(config['src_path'], w.animal, w.session, 'probe.json')):
        probe_path = os.path.join(config['src_path'], w.animal, w.session, 'probe.json')
    elif os.path.isfile(config['kilosort']['probe_path']):
        probe_path = config['kilosort']['probe_path']
    else:
        raise FileNotFoundError('No probe file found - cant detect electrodes')
    
    with open(probe_path, 'r') as f:
        probe = json.loads(f.read())

    electrodes = list(np.unique(np.array(probe['kcoords'], dtype=np.int16)) + 1)
    return expand(k_path(w.animal, w.session, 'NS.clu.' + '{electrode}'), electrode=electrodes)


rule do_kilosort:
    input:
        settings=ancient(os.path.join(config['src_path'], '{animal}', '{session}', 'kilosort.json')),
        probe=ancient(os.path.join(config['src_path'], '{animal}', '{session}', 'probe.json')),
        dat_file=ancient(k_path('{animal}', '{session}', '{session}.dat'))
    output:
        # put the whitened filtered path here
        k_path('{animal}', '{session}', 'spike_times.npy'),
        k_path('{animal}', '{session}', 'spike_clusters.npy'),
        k_path('{animal}', '{session}', 'templates.npy')
    conda:
        "/mnt/nevermind.data-share/ag-grothe/AG_Pecka/envs/kilosort"
    script:
        "../scripts/kilosort.py"


rule kilosort2neurosuite:
    input:
        probe=ancient(os.path.join(config['src_path'], '{animal}', '{session}', 'probe.json')),
        st=k_path('{animal}', '{session}', 'spike_times.npy'),
        sc=k_path('{animal}', '{session}', 'spike_clusters.npy'),
        tp=k_path('{animal}', '{session}', 'templates.npy')
    output:
        k_path('{animal}', '{session}', 'NS.clu.{electrode}')
    script:
        "../scripts/kilosort2neurosuite.py"


 # finalize processing
rule kilosort_ready:
    input:
        clu=get_clu_input_names
    output:
        k_path('{animal}', '{session}', 'kilosort.ready')
    params:
        session="{session}",
        animal="{animal}"
    shell:
        "touch %s" % k_path('{params.animal}', '{params.session}', 'kilosort.ready')