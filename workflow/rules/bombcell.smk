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
        templates=k_path("{animal}", "{session}", "{stream}", "templates.npy"),
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
        tmpdir="/tmp",
        bombcell=int(config.get("bombcell", {}).get("resource_slots", 1)),
        io_heavy=int(config.get("bombcell", {}).get("io_heavy", 1))
    script:
        "../scripts/bombcell.py"
