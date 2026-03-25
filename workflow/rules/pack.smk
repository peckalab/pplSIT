import os


def sync_none_path(wc):
    return os.path.join(config["dst_path"], wc.animal, wc.session, "sync", "sounds_sync.none.h5")

def sync_ephys_path(wc):
    return os.path.join(config["dst_path"], wc.animal, wc.session, "sync", "sounds_sync.ephys.h5")

# def wants_ephys_sync(wc) -> bool:
#     """
#     Decide whether the session *requests* ephys-based sync based on manual.json.
#     Based on your current convention:
#       - manual['ephys']['offset'] is int => manual only (no ephys sync)
#       - manual['ephys']['offset'] is dict => ephys sync requested (OneBox_ADC etc)
#     """
#     p = manual_json_path(wc)
#     with open(p) as f:
#         m = json.load(f)
#     off = m.get("ephys", {}).get("offset", 0)
#     return isinstance(off, dict)

# def chosen_sync_path(wc):
#     """
#     Choose ephys sync file if requested and ephys is staged,
#     else use the none placeholder.
#     """
#     if wants_ephys_sync(wc) and os.path.exists(ephys_staged_marker(wc)):
#         return sync_ephys_path(wc)
#     return sync_none_path(wc)

def has_raw_ephys_session(wc) -> bool:
    session_path = os.path.join(config["src_path"], wc.animal, wc.session)

    found_dat = False
    found_settings = False

    for dirpath, dirnames, filenames in os.walk(session_path):
        if "ephys" in dirnames:
            dirnames.remove("ephys")

        if "settings.xml" in filenames:
            found_settings = True
        if any(f.endswith(".dat") for f in filenames):
            found_dat = True

        if found_dat and found_settings:
            return True

    return False


def optional_file(path):
    return path if os.path.exists(path) else []


rule pack_base:
    input:
        positions=os.path.join(config["src_path"], "{animal}", "{session}", "positions.csv"),
        events=os.path.join(config["src_path"], "{animal}", "{session}", "events.csv"),
        sounds=os.path.join(config["src_path"], "{animal}", "{session}", "sounds.csv"),
        cfg=os.path.join(config["src_path"], "{animal}", "{session}", "{session}.json"),
        manual=os.path.join(config["src_path"], "{animal}", "{session}", "manual.json"),
        islands=lambda wc: optional_file(
            os.path.join(config["src_path"], wc.animal, wc.session, "islands.csv")
        )
    output:
        base=os.path.join(config["dst_path"], "{animal}", "{session}", "meta.base.h5")
    script:
        "../scripts/pack/base.py"


rule sounds_sync_none:
    output:
        sync=os.path.join(config["dst_path"], "{animal}", "{session}", "sync", "sounds_sync.none.h5")
    run:
        import os, h5py
        os.makedirs(os.path.dirname(output.sync), exist_ok=True)
        with h5py.File(output.sync, "w") as f:
            f.attrs["mode"] = "no_sync"


rule sounds_sync_ephys:
    input:
        staged=os.path.join(config["src_path"], "{animal}", "{session}", "ephys", ".STAGED"),
        settings=os.path.join(config["src_path"], "{animal}", "{session}", "ephys", "settings.xml"),
        manual=os.path.join(config["src_path"], "{animal}", "{session}", "manual.json"),
        sounds=os.path.join(config["src_path"], "{animal}", "{session}", "sounds.csv"),
        events=os.path.join(config["src_path"], "{animal}", "{session}", "events.csv"),
    output:
        sync=os.path.join(config["dst_path"], "{animal}", "{session}", "sync", "sounds_sync.ephys.h5")
    script:
        "../scripts/pack/sync_sounds.py"


rule pack_merge:
    input:
        base=os.path.join(config["dst_path"], "{animal}", "{session}", "meta.base.h5"),
        manual=os.path.join(config["src_path"], "{animal}", "{session}", "manual.json"),
        sync_none=os.path.join(config["dst_path"], "{animal}", "{session}", "sync", "sounds_sync.none.h5"),

        # trigger ephys sync only for sessions that have raw ephys data
        sync_ephys=lambda wc: (
            [os.path.join(config["dst_path"], wc.animal, wc.session, "sync", "sounds_sync.ephys.h5")]
            if has_raw_ephys_session(wc) else []
        ),
    output:
        meta=os.path.join(config["dst_path"], "{animal}", "{session}", "meta.h5")
    script:
        "../scripts/pack/merge.py"

# def conditional_input_isl(wildcards):
#     f = os.path.join(config['src_path'], wildcards.animal, wildcards.session, 'islands.csv')
#     return f if os.path.exists(f) else []

# # this function checks if continuous.dat exists at all, 
# # if so the final path of dat in the raw folder is used as input,
# # which will trigger its creation in copy.smk (if needed) as well as the adc files if it makes sense
# def sync_inputs(wildcards):
#     base = os.path.join(config['src_path'], wildcards.animal, wildcards.session)

#     for dirpath, _, filenames in os.walk(base):
#         if "continuous.dat" in filenames:
#             return [
#                 os.path.join(config['src_path'], '{animal}', '{session}', '{session}.dat'),
#                 os.path.join(config['src_path'], '{animal}', '{session}', 'timestamps.npy'),
#             ]

#     return []  # important: return [] if not found

# rule pack:
#     input:
#         positions = os.path.join(config['src_path'], '{animal}', '{session}', 'positions.csv'),
#         events = os.path.join(config['src_path'], '{animal}', '{session}', 'events.csv'),
#         sounds = os.path.join(config['src_path'], '{animal}', '{session}', 'sounds.csv'),
#         islands = conditional_input_isl,
#         cfg = os.path.join(config['src_path'], '{animal}', '{session}', '{session}' + '.json'),
#         manual = os.path.join(config['src_path'], '{animal}', '{session}', 'manual.json'),
        
#         optional_sync_inputs = sync_inputs,
#     output:
#         meta = os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5')
#     script:
#         "../scripts/pack.py"