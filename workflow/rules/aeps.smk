import os


rule aeps_lfp_metrics_stream:
    input:
        meta      = os.path.join(config["dst_path"], "{animal}", "{session}", "meta.h5"),
        lfp       = os.path.join(config["dst_path"], "{animal}", "{session}", "LFP", "{stream}", "lfp.h5"),
        baseline  = os.path.join(config["dst_path"], "{animal}", "{session}", "LFP", "{stream}", "baseline.h5"),
        artifacts = os.path.join(config["dst_path"], "{animal}", "{session}", "LFP", "{stream}", "artifacts.h5")
    output:
        aeps_lfp_metrics = os.path.join(config["dst_path"], "{animal}", "{session}", "AEP", "{stream}", "aeps_lfp_metrics.h5")
    script:
        "../scripts/aeps/aeps_lfp_metrics.py"


rule aeps_ITPC_stream:
    input:
        aeps = os.path.join(config["dst_path"], "{animal}", "{session}", "AEP", "{stream}", "AEPs.h5"),
        segm = os.path.join(config["dst_path"], "{animal}", "{session}", "analysis", "segments.h5")
    output:
        ITPC = os.path.join(config["dst_path"], "{animal}", "{session}", "AEP", "{stream}", "aeps_ITPC.h5")
    script:
        "../scripts/aeps/aeps_ITPC.py"


rule aeps_ITPC_plot_stream:
    input:
        ITPC = os.path.join(config["dst_path"], "{animal}", "{session}", "AEP", "{stream}", "aeps_ITPC.h5")
    output:
        ITPC_plots = os.path.join(config["dst_path"], "{animal}", "{session}", "AEP", "{stream}", "aeps_ITPC.pdf")
    script:
        "../scripts/aeps/aeps_ITPC_plot.py"


rule aeps_EV_SU_metrics_stream:
    input:
        aeps = os.path.join(config["dst_path"], "{animal}", "{session}", "AEP", "{stream}", "AEPs.h5")
    output:
        metrics = os.path.join(config["dst_path"], "{animal}", "{session}", "AEP", "{stream}", "aeps_EV_SU_metrics.h5")
    script:
        "../scripts/aeps/aeps_EV_SU_metrics.py"


rule plot_AEP_profiles_stream:
    input:
        meta = os.path.join(config["dst_path"], "{animal}", "{session}", "meta.h5"),
        aeps = os.path.join(config["dst_path"], "{animal}", "{session}", "AEP", "{stream}", "AEPs.h5")
    output:
        os.path.join(config["dst_path"], "{animal}", "{session}", "AEP", "{stream}", "aeps_profiles.pdf")
    script:
        "../scripts/aeps/aeps_profiles.py"


rule extract_aeps_stream:
    input:
        manual       = os.path.join(config["src_path"], "{animal}", "{session}", "manual.json"),
        meta         = os.path.join(config["dst_path"], "{animal}", "{session}", "meta.h5"),
        lfp          = os.path.join(config["dst_path"], "{animal}", "{session}", "LFP", "{stream}", "lfp.h5"),
        lfp_base     = os.path.join(config["dst_path"], "{animal}", "{session}", "LFP", "{stream}", "baseline.h5"),
        artifacts    = os.path.join(config["dst_path"], "{animal}", "{session}", "LFP", "{stream}", "artifacts.h5"),
        aeps_metrics = os.path.join(config["dst_path"], "{animal}", "{session}", "AEP", "{stream}", "aeps_lfp_metrics.h5"),
        init         = os.path.join(config["src_path"], "{animal}", "{session}", ".templates_initialized"),
    output:
        aeps = os.path.join(config["dst_path"], "{animal}", "{session}", "AEP", "{stream}", "AEPs.h5")
    script:
        "../scripts/aeps/aeps_extract.py"


rule aep_ready_session:
    input:
        staged = ephys_staged_marker,
        profiles = lambda wc: expand(
            os.path.join(config["dst_path"], wc.animal, wc.session, "AEP", "{stream}", "aeps_profiles.pdf"),
            stream=streams_for_session_wc(wc)
        ),
        ev_su_metrics = lambda wc: expand(
            os.path.join(config["dst_path"], wc.animal, wc.session, "AEP", "{stream}", "aeps_EV_SU_metrics.h5"),
            stream=streams_for_session_wc(wc)
        ),
        itpc_plots = lambda wc: expand(
            os.path.join(config["dst_path"], wc.animal, wc.session, "AEP", "{stream}", "aeps_ITPC.pdf"),
            stream=streams_for_session_wc(wc)
        )
    output:
        os.path.join(config["dst_path"], "{animal}", "{session}", "AEP", "aep.ready")
    shell:
        "touch {output}"
