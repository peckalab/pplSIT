import os, sys

def get_session_src_dir(config, animal, session):
    return os.path.join(config["src_path"], animal, session)

def get_session_dst_dir(config, animal, session):
    return os.path.join(config["dst_path"], animal, session)

def k_path(animal, session, *parts):
    """processed root -> dst/<animal>/<session>/kilosort[/...]"""
    return os.path.join(config["dst_path"], animal, session, "kilosort", *parts)

def lfp_path(animal, session, *parts):
    return os.path.join(config["dst_path"], animal, session, "LFP", *parts)

def manual_json_path(wc):
    return os.path.join(config["src_path"], wc.animal, wc.session, "manual.json")
    
def ephys_stream_dir(animal, session, stream):
    # your step-1 staging output
    return os.path.join(config["src_path"], animal, session, "ephys", stream)

def ephys_stream_dat(wc):
    """Find the single .dat file in src/<animal>/<session>/ephys/<stream>/."""
    d = os.path.join(config["src_path"], wc.animal, wc.session, "ephys", wc.stream)
    dats = sorted(glob.glob(os.path.join(d, "*.dat")))
    if len(dats) != 1:
        raise ValueError(f"Expected exactly one .dat in {d}, found: {dats}")
    return dats[0]

def ephys_staged_marker(wc):
    return os.path.join(config["src_path"], wc.animal, wc.session, "ephys", ".STAGED")

def kilo_staged_marker(wc):
    return os.path.join(config["dst_path"], wc.animal, wc.session, "kilosort", ".STAGED")

def streams_for_session_from_ephys(wc):
    raw_ephys_dir = os.path.dirname(ephys_staged_marker(wc))
    if not os.path.exists(raw_ephys_dir):
        return []
        
    streams = []
    for name in os.listdir(raw_ephys_dir):
        p = os.path.join(raw_ephys_dir, name)
        if not os.path.isdir(p):
            continue
        # keep only real probe streams
        if name.startswith("Probe"):
            streams.append(name)
    return sorted(streams)

def streams_for_session_from_kilo(wc):
    staged = os.path.join(config["dst_path"], wc.animal, wc.session, "kilosort", ".STAGED")
    if not os.path.exists(staged):
        return []
    with open(staged) as f:
        return [ln.strip() for ln in f if ln.strip()]


def is_real_ephys_stream(stream_name: str) -> bool:
    """
    Keep real neural streams, exclude obvious ADC streams.
    """
    return "adc" not in stream_name.lower()


def streams_for_session(config, animal, session):
    """
    Return stream names for a session.

    Prefer staged ephys/<stream>/ folders if present.
    Otherwise infer from the raw session tree using the same naming logic
    as stage_ephys_streams().
    """
    sess_dir = get_session_src_dir(config, animal, session)
    ephys_dir = os.path.join(sess_dir, "ephys")

    # Case 1: staged ephys exists
    if os.path.isdir(ephys_dir):
        staged_streams = []
        for name in os.listdir(ephys_dir):
            p = os.path.join(ephys_dir, name)
            if not os.path.isdir(p):
                continue
            if is_real_ephys_stream(name):
                staged_streams.append(name)

        if staged_streams:
            return sorted(staged_streams)

    # Case 2: infer from raw tree
    if not os.path.isdir(sess_dir):
        return []

    raw_streams = set()

    for dirpath, dirnames, filenames in os.walk(sess_dir):
        if "ephys" in dirnames:
            dirnames.remove("ephys")

        if os.path.abspath(dirpath) == os.path.abspath(sess_dir):
            continue

        dat_files = [f for f in filenames if f.endswith(".dat")]
        if not dat_files:
            continue

        parent_dirname = os.path.basename(dirpath)
        stream_name = stream_name_from_dir(parent_dirname)

        if is_real_ephys_stream(stream_name):
            raw_streams.add(stream_name)

    return sorted(raw_streams)


def streams_for_session_wc(wc):
    return streams_for_session(config, wc.animal, wc.session)