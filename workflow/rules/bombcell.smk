import os


def bombcell_ready_files(wc):
    streams = configured_streams_for_session_wc("bombcell", wc)
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
        raw_waveforms_ready=k_path("{animal}", "{session}", "{stream}", "bombcell", "RawWaveforms", "RawWaveforms.ready"),
        unit_type=k_path("{animal}", "{session}", "{stream}", "cluster_bc_unitType.tsv"),
        metrics_csv=k_path("{animal}", "{session}", "{stream}", "bombcell", "templates._bc_qMetrics.csv"),
        metrics_parquet=k_path("{animal}", "{session}", "{stream}", "bombcell", "templates._bc_qMetrics.parquet"),
    params:
        param_overrides=lambda wc: config.get("bombcell", {}).get("param_overrides", {}),
        unit_match_waveforms=lambda wc: config.get("bombcell", {}).get("unit_match_waveforms", False),
        force_unit_match_waveform_reextract=lambda wc: config.get("bombcell", {}).get("force_unit_match_waveform_reextract", False),
        unit_match_spike_width=lambda wc: config.get("unit_match", {}).get("spike_width", None),
        unit_match_sample_amount=lambda wc: config.get("unit_match", {}).get("sample_amount", None),
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
