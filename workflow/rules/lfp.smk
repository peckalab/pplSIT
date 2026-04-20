import os
import numpy as np


rule extract_lfp_raw_stream:
    input:
        staged=ephys_staged_marker
    output:
        os.path.join(config["dst_path"], "{animal}", "{session}", "LFP", "{stream}", "lfp.h5")
    script:
        "../scripts/lfp/lfp.py"


rule extract_lfp_artifacts_stream:
    input:
        ancient(os.path.join(config["dst_path"], "{animal}", "{session}", "LFP", "{stream}", "lfp.h5"))
    output:
        artifacts=os.path.join(config["dst_path"], "{animal}", "{session}", "LFP", "{stream}", "artifacts.h5"),
        artifacts_pdf=os.path.join(config["dst_path"], "{animal}", "{session}", "LFP", "{stream}", "artifacts.pdf")
    script:
        "../scripts/lfp/artifacts.py"


rule extract_lfp_baseline_stream:
    input:
        meta=os.path.join(config["dst_path"], "{animal}", "{session}", "meta.h5"),
        lfp_h5=ancient(os.path.join(config["dst_path"], "{animal}", "{session}", "LFP", "{stream}", "lfp.h5")),
        artifacts=os.path.join(config["dst_path"], "{animal}", "{session}", "LFP", "{stream}", "artifacts.h5")
    output:
        lfp_base=os.path.join(config["dst_path"], "{animal}", "{session}", "LFP", "{stream}", "baseline.h5"),
        lfp_base_plot=os.path.join(config["dst_path"], "{animal}", "{session}", "LFP", "{stream}", "baseline.pdf")
    script:
        "../scripts/lfp/baseline.py"


rule lfp_ready_session:
    input:
        staged=ephys_staged_marker,
        baselines=lambda wc: expand(
            os.path.join(config["dst_path"], wc.animal, wc.session, "LFP", "{stream}", "baseline.h5"),
            stream=streams_for_session_wc(wc)
        )
    output:
        os.path.join(config["dst_path"], "{animal}", "{session}", "LFP", "lfp.ready")
    shell:
        "touch {output}"