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
    params:
        session = "{session}",
        animal = "{animal}"
    script:
        "../scripts/units.py"


rule unit_response_metrics:
    input:
        meta =os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        units=os.path.join(config['dst_path'], '{animal}', '{session}', 'units.h5'),
        segm =os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'segments.h5')
    output:
        os.path.join(config['dst_path'], '{animal}', '{session}', 'units', 'unit_response_metrics.h5')
    script:
        "../scripts/units/unit_response_metrics.py"