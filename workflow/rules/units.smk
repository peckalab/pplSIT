import os

def guess_sorter(w):
    if config['units']['source'] == 'neurosuite':
        return n_path(w.animal, w.session, 'neurosuite.ready')
    else:
        return k_path(w.animal, w.session, 'kilosort.ready')

 
# dump units to HDF5
rule dump_units:
    input:
        meta=os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        clu=guess_sorter
    output:
        os.path.join(config['dst_path'], '{animal}', '{session}', 'units.h5')
    script:
        "../scripts/units.py"