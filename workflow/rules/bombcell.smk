import os


def bombcell_ready_files(wc):
    streams = streams_for_session_wc(wc)
    return expand(
        k_path(wc.animal, wc.session, "{stream}", "bombcell", "bombcell.ready"),
        stream=streams
    )


rule bombcell_stream:
    input:
        ready=k_path("{animal}", "{session}", "{stream}", "kilosort.ready"),
        st=k_path("{animal}", "{session}", "{stream}", "spike_times.npy"),
        sc=k_path("{animal}", "{session}", "{stream}", "spike_clusters.npy"),
        spike_templates=k_path("{animal}", "{session}", "{stream}", "spike_templates.npy"),
        templates=k_path("{animal}", "{session}", "{stream}", "templates.npy"),
        amplitudes=k_path("{animal}", "{session}", "{stream}", "amplitudes.npy"),
        whitening_mat_inv=k_path("{animal}", "{session}", "{stream}", "whitening_mat_inv.npy"),
        channel_positions=k_path("{animal}", "{session}", "{stream}", "channel_positions.npy"),
        dat=lambda wc: stream_dat_path(wc.animal, wc.session, wc.stream),
        settings=k_path("{animal}", "{session}", "{stream}", "settings.json"),
    output:
        ready=k_path("{animal}", "{session}", "{stream}", "bombcell", "bombcell.ready"),
        unit_type=k_path("{animal}", "{session}", "{stream}", "cluster_bc_unitType.tsv"),
        metrics_csv=k_path("{animal}", "{session}", "{stream}", "bombcell", "templates._bc_qMetrics.csv"),
        metrics_parquet=k_path("{animal}", "{session}", "{stream}", "bombcell", "templates._bc_qMetrics.parquet"),
    params:
        param_overrides=lambda wc: config.get("bombcell", {}).get("param_overrides", {}),
    threads: int(config.get("bombcell", {}).get("threads", 16))
    conda:
        config.get("bombcell", {}).get(
            "conda_env",
            "/mnt/nevermind.data-share/ag-grothe/AG_Pecka/envs/bombcell"
        )
    resources:
        tmpdir="/tmp"
    script:
        "../scripts/bombcell.py"
