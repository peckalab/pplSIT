import os, json, glob, sys
import numpy as np
import h5py
import xml.etree.ElementTree as ET

# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.sync import get_sound_events_from_ADC, get_sound_events_from_openephys


def _find_onebox(root: ET.Element) -> ET.Element:
    sc = root.find("SIGNALCHAIN")
    if sc is None:
        raise ValueError("SIGNALCHAIN not found in Open Ephys settings.xml")

    for proc in sc.findall("PROCESSOR"):
        if proc.attrib.get("name") == "OneBox" or proc.attrib.get("pluginName") == "OneBox":
            return proc
    raise ValueError("OneBox processor not found in Open Ephys settings.xml")


def _stream_info_from_settings(settings_xml: str, stream_name: str):
    """
    Return (sample_rate, channel_count) from OneBox STREAM attributes.
    """
    root = ET.parse(settings_xml).getroot()
    onebox = _find_onebox(root)

    for s in onebox.findall("STREAM"):
        if s.attrib.get("name") == stream_name:
            sr = float(s.attrib.get("sample_rate"))
            cc = int(s.attrib.get("channel_count"))
            return sr, cc

    raise ValueError(f"STREAM '{stream_name}' not found under OneBox in settings.xml")


def _single_file(glob_pattern: str, what: str) -> str:
    hits = sorted(glob.glob(glob_pattern))
    if len(hits) != 1:
        raise ValueError(f"Expected exactly one {what} matching {glob_pattern}, found: {hits}")
    return hits[0]


def _pick_reference_ephys_timestamps(ephys_root: str) -> str:
    """
    For OneBox_ADC sync we need ephys_ts_file to define time zero.
    Choose timestamps.npy from the first non-ADC stream (ProbeA/ProbeB).
    """
    if not os.path.isdir(ephys_root):
        raise ValueError(f"ephys_root does not exist: {ephys_root}")

    stream_dirs = sorted(
        d for d in os.listdir(ephys_root)
        if os.path.isdir(os.path.join(ephys_root, d))
    )

    probe_dirs = [d for d in stream_dirs if "adc" not in d.lower()]
    if not probe_dirs:
        raise ValueError(f"No non-ADC probe streams found under {ephys_root}")

    ts = os.path.join(ephys_root, probe_dirs[0], "timestamps.npy")
    if not os.path.exists(ts):
        raise FileNotFoundError(f"Reference ephys timestamps.npy not found: {ts}")

    return ts


def write_sync_h5(out_path: str, mode: str, manual: dict,
                  ev_detected: np.ndarray | None,
                  ev_synced: np.ndarray | None,
                  extra_attrs: dict | None = None):
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with h5py.File(out_path, "w") as f:
        f.attrs["mode"] = mode
        f.attrs["manual_offset_type"] = str(manual.get("ephys", {}).get("offset", {}).get("type", "unknown"))

        if extra_attrs:
            for k, v in extra_attrs.items():
                try:
                    f.attrs[k] = v
                except TypeError:
                    f.attrs[k] = str(v)

        if ev_detected is not None:
            f.create_dataset("events_detected", data=np.asarray(ev_detected))
        if ev_synced is not None:
            f.create_dataset("events_synced", data=np.asarray(ev_synced))


# ----------------- main -----------------

manual_path   = snakemake.input["manual"]
settings_xml  = snakemake.input["settings"]
sounds_csv    = snakemake.input["sounds"]
events_csv    = snakemake.input["events"]
ephys_root    = snakemake.input["ephys_root"]
out_sync_h5   = snakemake.output["sync"]

with open(manual_path, "r") as f:
    manual = json.load(f)

offset = manual.get("ephys", {}).get("offset", 0)

# If manual offset is not a dict => ephys sync not requested (this script should not be scheduled then)
if not isinstance(offset, dict):
    raise ValueError(
        f"manual['ephys']['offset'] is not a dict (got {type(offset)}). "
        "This sync_sounds.py is intended for dict-based ephys sync (e.g. OneBox_ADC)."
    )

sync_type = offset.get("type", None)
if sync_type is None:
    raise ValueError("manual['ephys']['offset']['type'] missing")

if sync_type == "OneBox_ADC":
    adc_channel = int(offset.get("channel", 11))
    event_th    = int(offset.get("threshold", 300))

    # ADC stream folder name produced by your staging: likely 'OneBox-ADC'
    adc_stream_name = "OneBox-ADC"
    adc_sr, adc_cc = _stream_info_from_settings(settings_xml, adc_stream_name)

    adc_dir = os.path.join(ephys_root, adc_stream_name)
    adc_dat = _single_file(os.path.join(adc_dir, "*.dat"), "ADC .dat")
    adc_ts  = os.path.join(adc_dir, "timestamps.npy")
    if not os.path.exists(adc_ts):
        raise FileNotFoundError(f"ADC timestamps.npy not found: {adc_ts}")

    ephys_ts = _pick_reference_ephys_timestamps(ephys_root)

    ev_detected, ev_synced = get_sound_events_from_ADC(
        adc_file=adc_dat,
        adc_ts_file=adc_ts,
        sounds_file=sounds_csv,
        events_file=events_csv,
        ephys_ts_file=ephys_ts,
        channel=adc_channel,
        event_th=event_th,
        s_rate=adc_sr,
        ch_no=adc_cc,
    )
    write_sync_h5(
        out_path=out_sync_h5,
        mode="OneBox_ADC",
        manual=manual,
        ev_detected=ev_detected,
        ev_synced=ev_synced,
        extra_attrs={
            "adc_channel": adc_channel,
            "threshold": event_th,
            "adc_sample_rate": float(adc_sr),
            "adc_channel_count": int(adc_cc),
            "adc_dat": adc_dat,
            "adc_ts": adc_ts,
            "ephys_ts": ephys_ts,
        }
    )

elif sync_type == "ephys":
    # NOTE: your current get_sound_events_from_openephys() in utils/sync.py expects a Neurosuite XML
    # (XMLHero), not Open Ephys settings.xml. If you still use legacy .xml in some sessions, you can
    # pass it via input and implement here. Otherwise fail loudly.
    raise ValueError(
        "sync_type 'ephys' currently expects a Neurosuite .xml (XMLHero) in utils/sync.py, "
        "but your pipeline now stages Open Ephys settings.xml. "
        "Either implement an Open Ephys version or switch to OneBox_ADC sync."
    )

else:
    raise ValueError(f"Unsupported sync type: {sync_type}")
