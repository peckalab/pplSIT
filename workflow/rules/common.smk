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
    p = ephys_staged_marker(wc)
    if not os.path.exists(p):
        return []
    with open(p) as f:
        return [ln.strip() for ln in f if ln.strip()]

def streams_for_session_from_kilo(wc):
    staged = os.path.join(config["dst_path"], wc.animal, wc.session, "kilosort", ".STAGED")
    if not os.path.exists(staged):
        return []
    with open(staged) as f:
        return [ln.strip() for ln in f if ln.strip()]
