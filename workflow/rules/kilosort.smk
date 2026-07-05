import os, sys, json
import numpy as np
import glob

from utils.probe import oe2kilosort


# ---------- helpers ----------

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


# ---------- the rules ----------

rule stage_streams_to_kilosort_folder:
    """
    Links dat + timestamps into:
      processed/<animal>/<session>/kilosort/<stream>/
    Creates per-stream probe.json.
    Writes marker file listing streams.
    """
    input:
        ephys_staged=ancient(os.path.join(config["src_path"], "{animal}", "{session}", "ephys", ".STAGED")),
        kilosort_json=ancient(config["kilosort"]["settings_path"]),
    output:
        marker=os.path.join(config["dst_path"], "{animal}", "{session}", "kilosort", ".STAGED"),
    run:
        import shutil
        animal = wildcards.animal
        session = wildcards.session

        streams = streams_for_session(config, animal, session)
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


if config.get("copy_kilosort_legacy", False):

    rule do_kilosort_stream:
        input:
            staged=ancient(lambda wc: staged_marker_path(wc.animal, wc.session))
        output:
            st=k_path("{animal}", "{session}", "{stream}", "spike_times.npy"),
            sc=k_path("{animal}", "{session}", "{stream}", "spike_clusters.npy"),
            tp=k_path("{animal}", "{session}", "{stream}", "templates.npy"),
        run:
            import os
            import shutil
            from pathlib import Path

            ks_root = Path(config["dst_path"]) / wildcards.animal / wildcards.session / "kilosort"
            stream_dir = ks_root / wildcards.stream
            stream_dir.mkdir(parents=True, exist_ok=True)

            exclude_names = {
                wildcards.stream,   # do not recurse into destination subfolder
                ".STAGED",
                "kilosort.ready",
            }

            def should_skip(path: Path) -> bool:
                name = path.name

                if name in exclude_names:
                    return True

                if path.is_file() and name.endswith(".dat"):
                    return True

                return False

            copied = []

            for item in ks_root.iterdir():
                if should_skip(item):
                    continue

                dst = stream_dir / item.name

                if item.is_file():
                    shutil.copy2(item, dst)
                    copied.append((str(item), str(dst)))

                elif item.is_dir():
                    if dst.exists():
                        shutil.rmtree(dst)
                    shutil.copytree(item, dst)
                    copied.append((str(item), str(dst)))

            # Ensure required declared outputs exist after copy
            missing = [p for p in output if not os.path.exists(p)]
            if missing:
                raise FileNotFoundError(
                    "Legacy Kilosort copy completed, but required outputs are missing: "
                    + ", ".join(missing)
                )

            print(f"Copied legacy Kilosort contents from {ks_root} -> {stream_dir}")
            for src, dst in copied:
                print(f"  {src} -> {dst}")

else:

    rule do_kilosort_stream:
        input:
            staged=ancient(lambda wc: staged_marker_path(wc.animal, wc.session))
        output:
            st=k_path("{animal}", "{session}", "{stream}", "spike_times.npy"),
            sc=k_path("{animal}", "{session}", "{stream}", "spike_clusters.npy"),
            tp=k_path("{animal}", "{session}", "{stream}", "templates.npy")
        conda:
            "/mnt/nevermind.data-share/ag-grothe/AG_Pecka/envs/kilosort"
        resources:
            kilosort_gpu=int(config.get("kilosort", {}).get("resource_slots", 1))
        script:
            "../scripts/kilosort.py"


rule kilosort_ready_stream:
    input:
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
            os.path.join(
                config["dst_path"],
                wc.animal,
                wc.session,
                "kilosort",
                "{stream}",
                "kilosort.ready"
            ),
            stream=streams_for_session_wc(wc),
        )
    output:
        ready=k_path("{animal}", "{session}", "kilosort.ready")
    shell:
        "touch {output}"
