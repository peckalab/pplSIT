import os


def read_streams_from_marker(marker_path):
    with open(marker_path, "r") as f:
        return [ln.strip() for ln in f if ln.strip()]

def ks_marker_from_checkpoint(wc):
    # This forces Snakemake to run the checkpoint first and then gives us the real output path.
    ckpt = checkpoints.stage_streams_to_kilosort_folder.get(animal=wc.animal, session=wc.session)
    return ckpt.output.marker  # this is .../kilosort/.STAGED

def ks_ready_files(wc):
    marker = ks_marker_from_checkpoint(wc)
    streams = read_streams_from_marker(marker)
    return expand(k_path(wc.animal, wc.session, "{stream}", "kilosort.ready"), stream=streams)

def guess_sorter_inputs(wc):
    if config["units"]["source"] == "neurosuite":
        return [n_path(wc.animal, wc.session, "neurosuite.ready")]
    else:
        marker = ks_marker_from_checkpoint(wc)
        return [marker] + ks_ready_files(wc)

rule dump_units:
    input:
        meta=os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        clu=guess_sorter_inputs
    output:
        os.path.join(config['dst_path'], '{animal}', '{session}', 'units.h5')
    params:
        session="{session}",
        animal="{animal}"
    script:
        "../scripts/units.py"


# def guess_sorter(w):
#     if config['units']['source'] == 'neurosuite':
#         return n_path(w.animal, w.session, 'neurosuite.ready')
#     else:
#         return k_path(w.animal, w.session, 'kilosort.ready')

 
# # dump units to HDF5
# rule dump_units:
#     input:
#         meta=os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
#         clu=guess_sorter
#     output:
#         os.path.join(config['dst_path'], '{animal}', '{session}', 'units.h5')
#     params:
#         session = "{session}",
#         animal = "{animal}"
#     script:
#         "../scripts/units.py"


rule unit_response_metrics:
    input:
        meta =os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        units=os.path.join(config['dst_path'], '{animal}', '{session}', 'units.h5'),
        segm =os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'segments.h5')
    output:
        os.path.join(config['dst_path'], '{animal}', '{session}', 'units', 'unit_response_metrics.h5')
    script:
        "../scripts/units/unit_response_metrics.py"