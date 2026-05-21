import os


rule aeps_lfp_metrics_stream:
    input:
        meta      = os.path.join(config["dst_path"], "{animal}", "{session}", "meta.h5"),
        lfp       = os.path.join(config["dst_path"], "{animal}", "{session}", "LFP", "{stream}", "lfp.h5"),
        baseline  = lfp_baseline_h5_for_wc,
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


rule plot_AEP_profiles_passive_stream:
    input:
        meta = os.path.join(config["dst_path"], "{animal}", "{session}", "meta.h5"),
        aeps = os.path.join(config["dst_path"], "{animal}", "{session}", "AEP", "{stream}", "AEPs.h5")
    output:
        os.path.join(config["dst_path"], "{animal}", "{session}", "AEP", "{stream}", "aeps_profiles_passive.pdf")
    script:
        "../scripts/aeps/aeps_profiles_passive.py"


rule extract_aeps_stream:
    input:
        meta         = os.path.join(config["dst_path"], "{animal}", "{session}", "meta.h5"),
        lfp          = os.path.join(config["dst_path"], "{animal}", "{session}", "LFP", "{stream}", "lfp.h5"),
        lfp_base     = lfp_baseline_h5_for_wc,
        artifacts    = os.path.join(config["dst_path"], "{animal}", "{session}", "LFP", "{stream}", "artifacts.h5"),
        aeps_metrics = os.path.join(config["dst_path"], "{animal}", "{session}", "AEP", "{stream}", "aeps_lfp_metrics.h5"),
        init         = os.path.join(config["src_path"], "{animal}", "{session}", ".templates_initialized"),
    params:
        manual=os.path.join(config["src_path"], "{animal}", "{session}", "manual.json"),
    output:
        aeps = os.path.join(config["dst_path"], "{animal}", "{session}", "AEP", "{stream}", "AEPs.h5")
    script:
        "../scripts/aeps/aeps_extract.py"


rule aep_ready_session:
    input:
        staged = ephys_staged_marker,
        profiles = aep_profile_paths_for_session_wc,
        ev_su_metrics = lambda wc: expand(
            os.path.join(config["dst_path"], wc.animal, wc.session, "AEP", "{stream}", "aeps_EV_SU_metrics.h5"),
            stream=streams_for_session_wc(wc)
        ),
        itpc_plots = aep_itpc_plot_paths_for_session_wc
    output:
        os.path.join(config["dst_path"], "{animal}", "{session}", "AEP", "aep.ready")
    shell:
        "touch {output}"
