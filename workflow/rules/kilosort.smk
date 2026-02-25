import os, sys, json
import numpy as np
import glob

from utils.probe import oe2kilosort


# ---------- helpers ----------

def list_streams_for_kilosort(config, animal, session):
    """
    Return stream names under src/<animal>/<session>/ephys/*, excluding ADC streams.
    """
    ephys_root = os.path.join(get_session_src_dir(config, animal, session), "ephys")
    if not os.path.isdir(ephys_root):
        return []

    streams = []
    for name in sorted(os.listdir(ephys_root)):
        full = os.path.join(ephys_root, name)
        if not os.path.isdir(full):
            continue
        if "adc" in name.lower():
            continue
        streams.append(name)
    return streams

def safe_hardlink(src: str, dst: str, overwrite: bool = True) -> None:
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    if overwrite and os.path.lexists(dst):
        os.remove(dst)
    os.link(src, dst)

def write_json_atomic(data: dict, path: str) -> None:
    """
    Write JSON atomically (tmp -> rename) to avoid partial files if interrupted.
    """
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(data, f, indent=2)
        f.write("\n")
    os.replace(tmp, path)

def get_settings_xml_path(config, animal, session) -> str:
    """
    settings.xml is staged by step-1 into src/<animal>/<session>/ephys/settings.xml
    """
    p = os.path.join(get_session_src_dir(config, animal, session), "ephys", "settings.xml")
    if not os.path.exists(p):
        raise ValueError(f"settings.xml not found at expected path: {p}")
    return p

def stage_stream_to_processed(config, animal, session, stream):
    """
    Hardlink dat + timestamps for one stream into
    dst/<animal>/<session>/kilosort/<stream>/...
    AND write probe.json using oe2kilosort(settings.xml, stream).
    """
    src_stream_dir = os.path.join(get_session_src_dir(config, animal, session), "ephys", stream)
    dst_stream_dir = os.path.join(get_session_dst_dir(config, animal, session), "kilosort", stream)

    # Find dat files inside src stream dir
    dats = sorted(glob.glob(os.path.join(src_stream_dir, "*.dat")))
    if len(dats) == 0:
        raise ValueError(f"No .dat found in {src_stream_dir}")
    if len(dats) > 1:
        raise ValueError(f"Multiple .dat files found in {src_stream_dir}: {dats}")

    # Link dat
    src_dat = dats[0]
    dst_dat = os.path.join(dst_stream_dir, os.path.basename(src_dat))
    safe_hardlink(src_dat, dst_dat, overwrite=True)

    # Link timestamps if present
    src_ts = os.path.join(src_stream_dir, "timestamps.npy")
    if os.path.exists(src_ts):
        dst_ts = os.path.join(dst_stream_dir, "timestamps.npy")
        safe_hardlink(src_ts, dst_ts, overwrite=True)

    # Write probe.json
    settings_xml = get_settings_xml_path(config, animal, session)
    probe_cfg = oe2kilosort(settings_xml, stream)  # <-- your updated oe2kilosort(xml, probe_name)
    probe_json_path = os.path.join(dst_stream_dir, "probe.json")
    write_json_atomic(probe_cfg, probe_json_path)

def staged_marker_path(animal, session):
    return k_path(animal, session, ".STAGED")

def stream_dat_path(animal, session, stream):
    """
    Return the single .dat inside dst/.../kilosort/<stream>/.
    This avoids hardcoding continuous.dat vs raw.dat.
    """
    stream_dir = k_path(animal, session, stream)
    dats = sorted(glob.glob(os.path.join(stream_dir, "*.dat")))
    if len(dats) != 1:
        raise ValueError(f"Expected exactly one .dat in {stream_dir}, found: {dats}")
    return dats[0]

def stream_probe_path(animal, session, stream):
    return k_path(animal, session, stream, "probe.json")

def session_kilosort_settings_path(animal, session, stream):
    return k_path(animal, session, stream, "settings.json")

def read_streams_from_marker(marker_path):
    if not os.path.exists(marker_path):
        return []
    with open(marker_path, "r") as f:
        return [ln.strip() for ln in f if ln.strip()]

def streams_after_staging(wc):
    """
    This is the key Fix A function.
    It forces Snakemake to run the checkpoint first.
    """
    ck = checkpoints.stage_streams_to_kilosort_folder.get(
        animal=wc.animal,
        session=wc.session,
    )
    marker = ck.output.marker
    
    streams = read_streams_from_marker(marker)

    if not streams:
        raise ValueError(f"No streams found in marker {marker}")

    return streams

# ---------- the rules ----------

checkpoint stage_streams_to_kilosort_folder:
    """
    Links dat + timestamps into:
      processed/<animal>/<session>/kilosort/<stream>/
    Creates per-stream probe.json.
    Writes marker file listing streams.
    """
    input:
        settings=lambda wc: os.path.join(config["src_path"], wc.animal, wc.session, "ephys", "settings.xml"),
        ephys_staged=os.path.join(config["src_path"], "{animal}", "{session}", "ephys", ".STAGED"),
        #kilosort_json=os.path.join(config["src_path"], "{animal}", "{session}", "kilosort.json")
        kilosort_json=ancient(config["kilosort"]["settings_path"]),
    output:
        marker=os.path.join(config["dst_path"], "{animal}", "{session}", "kilosort", ".STAGED"),
    run:
        import shutil
        animal = wildcards.animal
        session = wildcards.session

        streams = list_streams_for_kilosort(config, animal, session)
        if len(streams) == 0:
            raise ValueError(
                f"No non-ADC streams found in the raw. "
                f"Expected subfolders like ephys/ProbeA, ephys/ProbeB."
            )

        for stream in streams:
            stage_stream_to_processed(config, animal, session, stream)

        # Ensure kilosort folder exists
        ks_folder = os.path.dirname(output.marker)
        os.makedirs(ks_folder, exist_ok=True)

        # Copy kilosort settings into each stream folder as settings.json
        for stream in streams:
            stream_folder = os.path.join(ks_folder, stream)  # .../kilosort/<stream>
            os.makedirs(stream_folder, exist_ok=True)
            shutil.copy2(input.kilosort_json, os.path.join(stream_folder, "settings.json"))

        with open(output.marker, "w") as f:
            f.write("\n".join(streams) + "\n")


rule do_kilosort_stream:
    input:
        staged=lambda wc: staged_marker_path(wc.animal, wc.session),   # ensures staging happened
        settings=lambda wc: session_kilosort_settings_path(wc.animal, wc.session, wc.stream),
        probe=lambda wc: stream_probe_path(wc.animal, wc.session, wc.stream),
        dat=lambda wc: stream_dat_path(wc.animal, wc.session, wc.stream),
    output:
        st=k_path("{animal}", "{session}", "{stream}", "spike_times.npy"),
        sc=k_path("{animal}", "{session}", "{stream}", "spike_clusters.npy"),
        tp=k_path("{animal}", "{session}", "{stream}", "templates.npy"),
    conda:
        "/mnt/nevermind.data-share/ag-grothe/AG_Pecka/envs/kilosort"
    script:
        "../scripts/kilosort.py"


rule kilosort_ready_stream:
    input:
        probe=k_path("{animal}", "{session}", "{stream}", "probe.json"),
        st=k_path("{animal}", "{session}", "{stream}", "spike_times.npy"),
        sc=k_path("{animal}", "{session}", "{stream}", "spike_clusters.npy"),
        tp=k_path("{animal}", "{session}", "{stream}", "templates.npy"),
    output:
        ready=k_path("{animal}", "{session}", "{stream}", "kilosort.ready")
    shell:
        "touch {output.ready}"


rule kilosort_ready_session:
    input:
        staged=lambda wc: os.path.join(
            config["dst_path"], wc.animal, wc.session, "kilosort", ".STAGED"
        ),
        ready_files=lambda wc: expand(
            os.path.join(config["dst_path"], wc.animal, wc.session, "kilosort", "{stream}", "kilosort.ready"),
            stream=streams_after_staging(wc),
        )
    output:
        ready=k_path("{animal}", "{session}", "kilosort.ready")
    shell:
        "touch {output}"


# rule do_kilosort:
#     input:
#         settings=ancient(os.path.join(config['src_path'], '{animal}', '{session}', 'kilosort.json')),
#         probe=ancient(os.path.join(config['src_path'], '{animal}', '{session}', 'probe.json')),
#         dat_file=ancient(k_path('{animal}', '{session}', '{session}.dat'))
#     output:
#         # put the whitened filtered path here
#         k_path('{animal}', '{session}', 'spike_times.npy'),
#         k_path('{animal}', '{session}', 'spike_clusters.npy'),
#         k_path('{animal}', '{session}', 'templates.npy')
#     conda:
#         "/mnt/nevermind.data-share/ag-grothe/AG_Pecka/envs/kilosort"
#     script:
#         "../scripts/kilosort.py"


#  # finalize processing
# rule kilosort_ready:
#     input:
#         probe=k_path('{animal}', '{session}', 'probe.json'),
#         st=k_path('{animal}', '{session}', 'spike_times.npy'),
#         sc=k_path('{animal}', '{session}', 'spike_clusters.npy'),
#         tp=k_path('{animal}', '{session}', 'templates.npy')
#     output:
#         k_path('{animal}', '{session}', 'kilosort.ready')
#     params:
#         session="{session}",
#         animal="{animal}"
#     shell:
#         "touch %s" % k_path('{params.animal}', '{params.session}', 'kilosort.ready')