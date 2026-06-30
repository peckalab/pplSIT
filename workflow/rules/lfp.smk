import os
import numpy as np


rule extract_lfp_raw_stream:
    input:
        staged=ephys_staged_marker
    output:
        os.path.join(config["dst_path"], "{animal}", "{session}", "LFP", "{stream}", "lfp.h5")
    threads:
        int(config.get("lfp", {}).get("scheduler_threads", 16))
    resources:
        io_heavy=int(config.get("lfp", {}).get("io_heavy", 1))
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


rule extract_lfp_baseline_passive_stream:
    input:
        meta=os.path.join(config["dst_path"], "{animal}", "{session}", "meta.h5"),
        lfp_h5=ancient(os.path.join(config["dst_path"], "{animal}", "{session}", "LFP", "{stream}", "lfp.h5")),
        artifacts=os.path.join(config["dst_path"], "{animal}", "{session}", "LFP", "{stream}", "artifacts.h5")
    output:
        lfp_base=os.path.join(config["dst_path"], "{animal}", "{session}", "LFP", "{stream}", "baseline_passive.h5"),
        lfp_base_plot=os.path.join(config["dst_path"], "{animal}", "{session}", "LFP", "{stream}", "baseline_passive.pdf")
    script:
        "../scripts/lfp/baseline_passive.py"


rule lfp_ready_session:
    input:
        staged=ephys_staged_marker,
        baselines=lfp_baseline_h5_paths_for_session_wc
    output:
        os.path.join(config["dst_path"], "{animal}", "{session}", "LFP", "lfp.ready")
    shell:
        "touch {output}"
