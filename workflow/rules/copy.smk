import os


def session_dir(wc):
    return os.path.join(config["src_path"], wc.animal, wc.session)

def stream_name_from_dir(parent_dirname: str) -> str:
    """Derive stream name from the folder name containing the .dat."""
    if "." in parent_dirname and 'probe' in parent_dirname.lower():
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
    Crawl session_path and for each *.dat file found:
      - derive stream_name from the folder containing the .dat
      - hardlink dat into session_path/ephys/<stream_name>/<original_dat_filename>
      - hardlink timestamps.npy (if present) into the same folder

    Additionally:
      - find the first settings.xml anywhere under session_path and hardlink it to
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
        raise ValueError(f"There should be at least one .dat file in: {session_path}")


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
        #xml_t=ancient(config["template_xml"]),
        man_t=ancient(config["template_manual_json"]),
        ks_t=ancient(config["kilosort"]["settings_path"]),
    output:
        #xml=os.path.join(config["src_path"], "{animal}", "{session}", "{session}.xml"),
        man=os.path.join(config["src_path"], "{animal}", "{session}", "manual.json"),
        ks=os.path.join(config["src_path"], "{animal}", "{session}", "kilosort.json"),
    shell:
        # mkdir -p {session_dir(wildcards)}
        r"""
        cp {input.man_t} {output.man}
        cp {input.ks_t}  {output.ks}
        """


rule copy_ephys_ns:
    input:
        xml=ancient(os.path.join(config['src_path'], '{animal}', '{session}', '{session}' + '.xml')),
        dat=ancient(os.path.join(config['src_path'], '{animal}', '{session}', '{session}' + '.dat'))
    output:
        xml=n_path('{animal}', '{session}', '{session}.xml'),
        dat=n_path('{animal}', '{session}', '{session}.dat')
    shell:
        "ln {input.xml} {output.xml}; ln {input.dat} {output.dat}"


# rule create_xml_from_template:
#     input:
#         template=ancient(config['template_xml'])
#     output:
#         xml=os.path.join(config['src_path'], '{animal}', '{session}', '{session}' + '.xml')
#     shell:
#         "cp {input.template} {output.xml}"


# rule create_manual_json_from_template:
#     input:
#         template=ancient(config['template_manual_json'])
#     output:
#         man_json=os.path.join(config['src_path'], '{animal}', '{session}', 'manual.json')
#     shell:
#         "cp {input.template} {output.man_json}"


# rule create_kilosort_settings_from_template:
#     input:
#         template=ancient(config['kilosort']['settings_path'])
#     output:
#         kilo=os.path.join(config['src_path'], '{animal}', '{session}', 'kilosort.json')
#     shell:
#         "cp {input.template} {output.kilo}"

# rule move_dat_from_subfolder:
#     output:
#         dat=os.path.join(config['src_path'], '{animal}', '{session}', '{session}.dat'),
#         ts_path = os.path.join(config['src_path'], '{animal}', '{session}', 'timestamps.npy')
#     run:
#         session_path = os.path.join(config['src_path'], wildcards.animal, wildcards.session)
#         stage_ephys_streams(session_path)


# rule copy_ephys_ks:
#     input:
#         dat=ancient(os.path.join(config['src_path'], '{animal}', '{session}', '{session}' + '.dat'))
#     output:
#         dat=k_path('{animal}', '{session}', '{session}.dat')
#     shell:
#         "ln {input.dat} {output.dat}"

# Creates a hard link in the session folder to the actual raw files with recorded ephys data.
# For Neuropixels recordings, also creates a hard link to the ADC raw data file.
# rule move_dat_from_subfolder:
#     output:
#         dat=os.path.join(config['src_path'], '{animal}', '{session}', '{session}.dat'),
#         ts_path = os.path.join(config['src_path'], '{animal}', '{session}', 'timestamps.npy')
#     run:
#         import os
#         import shutil

#         # Define the source path
#         session_path = os.path.join(config['src_path'], wildcards.animal, wildcards.session)

#         dat_path = None
#         for dirpath, dirnames, filenames in os.walk(session_path):
#             for filename in [f for f in filenames if f.endswith('.dat')]:
#                 parent_dirname = os.path.basename(dirpath)
#                 dat_path = os.path.join(dirpath, filename)

#                 if parent_dirname.find('OneBox-ADC') > 0:  # this is ADC dat file, special case for NP
#                     subprocess.run(['ln', dat_path, os.path.join(session_path, 'ADC.dat')])

#                     # assume here should be timestamps file too - need to move it up as well
#                     adc_ts_path = os.path.join(dirpath, 'timestamps.npy')
#                     subprocess.run(['ln', adc_ts_path, os.path.join(session_path, 'ADC_timestamps.npy')])
#                 else:
#                     subprocess.run(['ln', dat_path, output.dat])
#                     ts_path = os.path.join(dirpath, 'timestamps.npy')
#                     subprocess.run(['ln', ts_path, os.path.join(session_path, 'timestamps.npy')])

#         if dat_path is None:
#             raise ValueError("There should be at least one .dat file in the session path")
