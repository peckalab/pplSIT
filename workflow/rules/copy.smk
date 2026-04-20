import os


def session_dir(wc):
    return os.path.join(config["src_path"], wc.animal, wc.session)

def stream_name_from_dir(parent_dirname: str) -> str:
    """Derive stream name from the folder name containing the .dat."""
    if "." in parent_dirname and ('probe' in parent_dirname.lower() or 'adc' in parent_dirname.lower()):
        return parent_dirname.split(".")[-1]
    return parent_dirname


def safe_hardlink(src: str, dst: str, overwrite: bool = True) -> None:
    """Create a hardlink dst -> src. Optionally overwrite dst if it exists."""
    os.makedirs(os.path.dirname(dst), exist_ok=True)

    if os.path.abspath(src) == os.path.abspath(dst):
        return

    if overwrite and os.path.lexists(dst):
        # lexists() also returns True for broken symlinks
        os.remove(dst)

    os.link(src, dst)


def stage_ephys_streams(session_path: str) -> None:
    """
    Crawl all subfolders under session_path (but not session_path itself) and for each *.dat file found:
      - derive stream_name from the folder containing the .dat
      - hardlink dat into session_path/ephys/<stream_name>/<original_dat_filename>
      - hardlink timestamps.npy (if present) into the same folder

    Additionally:
      - find the first settings.xml anywhere under session_path subfolders and hardlink it to
        session_path/ephys/settings.xml
    """
    found_any_dat = False
    staged_settings = False

    ephys_root = os.path.join(session_path, "ephys")
    os.makedirs(ephys_root, exist_ok=True)

    for dirpath, dirnames, filenames in os.walk(session_path):
        # Prune destination folder so we don't re-stage what we just created
        if "ephys" in dirnames:
            dirnames.remove("ephys")

        # Skip files directly in session_path itself; only process subfolders
        if os.path.abspath(dirpath) == os.path.abspath(session_path):
            continue

        # (1) Stage the first settings.xml we encounter
        if (not staged_settings) and ("settings.xml" in filenames):
            src_xml = os.path.join(dirpath, "settings.xml")
            dst_xml = os.path.join(ephys_root, "settings.xml")
            safe_hardlink(src_xml, dst_xml, overwrite=True)
            staged_settings = True

        # (2) Stage dat + timestamps
        dat_files = [f for f in filenames if f.endswith(".dat")]
        if not dat_files:
            continue

        parent_dirname = os.path.basename(dirpath)
        stream_name = stream_name_from_dir(parent_dirname)

        for dat_file in dat_files:
            found_any_dat = True

            src_dat = os.path.join(dirpath, dat_file)
            dst_dir = os.path.join(ephys_root, stream_name)
            dst_dat = os.path.join(dst_dir, dat_file)

            safe_hardlink(src_dat, dst_dat, overwrite=True)

            src_ts = os.path.join(dirpath, "timestamps.npy")
            if os.path.exists(src_ts):
                dst_ts = os.path.join(dst_dir, "timestamps.npy")
                safe_hardlink(src_ts, dst_ts, overwrite=True)

    if not found_any_dat:
        raise ValueError(f"There should be at least one .dat file in subfolders of: {session_path}")


rule stage_ephys_from_openephys_tree:
    """
    Creates hardlinks for all discovered streams into:
      src/<animal>/<session>/ephys/<stream>/{*.dat,timestamps.npy}
    Also stages the first settings.xml into:
      src/<animal>/<session>/ephys/settings.xml
    """
    output:
        marker=os.path.join(config["src_path"], "{animal}", "{session}", "ephys", ".STAGED"),
        ephys_xml=os.path.join(config["src_path"], "{animal}", "{session}", "ephys", "settings.xml")
    run:
        session_path = os.path.join(config["src_path"], wildcards.animal, wildcards.session)
        stage_ephys_streams(session_path)

        # write marker
        os.makedirs(os.path.dirname(output.marker), exist_ok=True)
        with open(output.marker, "w") as f:
            f.write("OK\n")


rule init_session_templates:
    input:
        man_t=ancient(config["template_manual_json"]),
        #ks_t=ancient(config["kilosort"]["settings_path"]),
    output:
        marker=os.path.join(config["src_path"], "{animal}", "{session}", ".templates_initialized")
    run:
        import os
        import shutil

        session_dir = os.path.dirname(output.marker)
        man = os.path.join(session_dir, "manual.json")
        #ks = os.path.join(session_dir, "kilosort.json")

        os.makedirs(session_dir, exist_ok=True)

        if not os.path.exists(man):
            shutil.copy2(input.man_t, man)

        # if not os.path.exists(ks):
        #     shutil.copy2(input.ks_t, ks)

        with open(output.marker, "w") as f:
            f.write("OK\n")


rule copy_ephys_ns:
    input:
        xml=ancient(os.path.join(config['src_path'], '{animal}', '{session}', '{session}' + '.xml')),
        dat=ancient(os.path.join(config['src_path'], '{animal}', '{session}', '{session}' + '.dat'))
    output:
        xml=n_path('{animal}', '{session}', '{session}.xml'),
        dat=n_path('{animal}', '{session}', '{session}.dat')
    shell:
        "ln {input.xml} {output.xml}; ln {input.dat} {output.dat}"
