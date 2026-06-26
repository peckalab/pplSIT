import os


def ks_marker_path(wc):
    return os.path.join(
        config["dst_path"], wc.animal, wc.session, "kilosort", ".STAGED"
    )


def ks_ready_files(wc):
    streams = streams_for_session_wc(wc)
    return expand(
        k_path(wc.animal, wc.session, "{stream}", "kilosort.ready"),
        stream=streams
    )


def guess_sorter_inputs(wc):
    if config["units"]["source"] == "neurosuite":
        return [n_path(wc.animal, wc.session, "neurosuite.ready")]
    else:
        inputs = [ks_marker_path(wc)] + ks_ready_files(wc)
        if config.get("units", {}).get("label_source", "ks") == "bombcell":
            inputs += bombcell_ready_files(wc)
        return inputs


rule dump_units:
    input:
        meta=os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        clu=guess_sorter_inputs,
        ephys_staged=os.path.join(config["src_path"], "{animal}", "{session}", "ephys", ".STAGED")
    output:
        os.path.join(config['dst_path'], '{animal}', '{session}', 'units.h5')
    params:
        session="{session}",
        animal="{animal}"
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
