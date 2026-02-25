import numpy as np
import os, json
import xml.etree.ElementTree as ET


def _find_onebox_processor(root: ET.Element) -> ET.Element:
    sc = root.find("SIGNALCHAIN")
    if sc is None:
        raise ValueError("SIGNALCHAIN not found in Open Ephys XML")

    # your legacy used dict(x.items())['name'] == 'OneBox'
    # keep that behavior (Open Ephys uses name + pluginName in different versions)
    for proc in sc.findall("PROCESSOR"):
        attrs = proc.attrib
        if attrs.get("name") == "OneBox" or attrs.get("pluginName") == "OneBox":
            return proc
    raise ValueError("OneBox processor not found in the OpenEphys signal chain")

def _dock_from_probe_name(probe_name: str) -> str:
    """
    Default mapping: ProbeA->'1', ProbeB->'2', ProbeC->'3', ...
    You can replace with an explicit dict if you ever rename streams.
    """
    if not probe_name.startswith("Probe") or len(probe_name) < 6:
        raise ValueError(f"Unexpected probe_name '{probe_name}'. Expected like 'ProbeA'/'ProbeB'.")
    letter = probe_name[-1].upper()
    dock = ord(letter) - ord("A") + 1
    if dock < 1:
        raise ValueError(f"Could not infer dock from probe_name '{probe_name}'")
    return str(dock)

def _get_channel_map_indices_for_stream(root: ET.Element, stream_name: str):
    """
    Try to pull channel map indices for the given stream from the 'Channel Map' processor.
    Falls back to the first stream if stream-specific selection isn't possible.
    Returns list[int].
    """
    sc = root.find("SIGNALCHAIN")
    if sc is None:
        raise ValueError("SIGNALCHAIN not found in Open Ephys XML")

    ch_map_proc = None
    for proc in sc.findall("PROCESSOR"):
        attrs = proc.attrib
        if attrs.get("name") == "Channel Map" or attrs.get("pluginName") == "Channel Map":
            ch_map_proc = proc
            break
    if ch_map_proc is None:
        raise ValueError("Channel Map processor not found in the OpenEphys signal chain")

    cust = ch_map_proc.find("CUSTOM_PARAMETERS")
    if cust is None or len(list(cust)) == 0:
        raise ValueError("Channel Map/CUSTOM_PARAMETERS not found or empty")

    # Structure in many Open Ephys XMLs:
    # CUSTOM_PARAMETERS
    #   <STREAM name="ProbeA"> <CHANNEL index="..."/> ... </STREAM>
    #   <STREAM name="ProbeB"> ...
    # In some, it's unnamed and you only get "first stream".
    streams = list(cust)

    # Prefer a stream element whose attrib name matches
    chosen = None
    for s in streams:
        if s.attrib.get("name") == stream_name:
            chosen = s
            break
    if chosen is None:
        # fallback: first stream-like element
        chosen = streams[0]

    indices = []
    for ch in list(chosen):
        if "index" in ch.attrib:
            indices.append(int(ch.attrib["index"]))
    return indices

def oe2kilosort(oe_xml_path: str, probe_name: str):
    """
    Extract probe information from Open Ephys settings XML and return a Kilosort-style dict.

    Args:
      oe_xml_path: path to settings.xml
      probe_name: e.g. 'ProbeA', 'ProbeB' (stream name)

    Behavior:
      - ProbeA -> NP_PROBE dock="1"
      - ProbeB -> NP_PROBE dock="2"
      - pulls ELECTRODE_XPOS/YPOS + CHANNELS (kcoords)
      - applies Channel Map reorder for the matching stream if available
      - drops channel index 384 if present (typical sync channel in 385-ch streams)
    """
    root = ET.parse(oe_xml_path).getroot()
    onebox = _find_onebox_processor(root)

    # sampling rate: pull from the matching stream in OneBox
    stream_elem = None
    for s in onebox.findall("STREAM"):
        if s.attrib.get("name") == probe_name:
            stream_elem = s
            break
    if stream_elem is None:
        raise ValueError(f"Stream '{probe_name}' not found under OneBox STREAM entries")

    sample_rate = float(stream_elem.attrib.get("sample_rate", "0") or 0)
    if sample_rate <= 0:
        # not fatal for chanMap, but usually indicates XML mismatch
        raise ValueError(f"Invalid/missing sample_rate for stream '{probe_name}'")

    editor = onebox.find("EDITOR")
    if editor is None:
        raise ValueError("OneBox/EDITOR section not found")

    dock = _dock_from_probe_name(probe_name)

    # pick NP_PROBE by dock
    probe_elem = None
    for p in editor.findall("NP_PROBE"):
        if p.attrib.get("dock") == dock:
            probe_elem = p
            break
    if probe_elem is None:
        available = [p.attrib.get("dock") for p in editor.findall("NP_PROBE")]
        raise ValueError(f"NP_PROBE with dock='{dock}' not found (available docks: {available})")

    channels = probe_elem.find("CHANNELS")
    x_pos = probe_elem.find("ELECTRODE_XPOS")
    y_pos = probe_elem.find("ELECTRODE_YPOS")
    if channels is None or x_pos is None or y_pos is None:
        raise ValueError(f"Missing CHANNELS/ELECTRODE_XPOS/ELECTRODE_YPOS for dock={dock}")

    # Ensure deterministic CH0..CH383 ordering (don't rely on dict .values() order)
    ch_keys = list(channels.attrib.keys())   # like ["CH0","CH1",...]
    if not ch_keys:
        raise ValueError(f"No CHANNELS entries found for dock={dock}")

    # Sort by numeric channel id
    ch_nums = sorted([int(k[2:]) for k in ch_keys])
    ch_count = len(ch_nums)

    # coords
    xc = [float(x_pos.attrib[f"CH{i}"]) for i in ch_nums]
    yc = [float(y_pos.attrib[f"CH{i}"]) for i in ch_nums]

    # kcoords: description like "...:X" in your legacy
    # keep same parse but make it robust
    kcoords = []
    for i in ch_nums:
        desc = channels.attrib.get(f"CH{i}", "")
        # old: desc.split(':')[1]
        if ":" in desc:
            kcoords.append(int(desc.split(":")[-1]))
        else:
            # fallback: single shank
            kcoords.append(0)

    # channel map indices from Channel Map processor
    channel_map_list = _get_channel_map_indices_for_stream(root, probe_name)

    # drop sync channel index 384 if present
    channel_map_list = [i for i in channel_map_list if i != 384]

    # Map channel_map_list into indices of ch_nums
    # channels_probe is typically [0..383]; but compute robustly:
    channels_probe = ch_nums
    channels_probe_index = {ch: idx for idx, ch in enumerate(channels_probe)}

    idxs_channels = []
    for ch in channel_map_list:
        if ch in channels_probe_index:
            idxs_channels.append(channels_probe_index[ch])
        # else ignore any weird indices (safety)

    if not idxs_channels:
        # If Channel Map doesn’t match, fall back to identity
        idxs_channels = list(range(ch_count))

    return {
        # for Kilosort this is usually 0..n-1 after reordering
        "chanMap": list(range(len(idxs_channels))),
        "xc": [xc[i] for i in idxs_channels],
        "yc": [yc[i] for i in idxs_channels],
        "kcoords": [kcoords[i] for i in idxs_channels],
        "n_chan": len(idxs_channels),
        "sample_rate": sample_rate,
        "probe_name": probe_name,
        "dock": dock,
        "probe_serial_number": probe_elem.attrib.get("probe_serial_number"),
        "headstage_serial_number": probe_elem.attrib.get("headstage_serial_number"),
    }


def kilosort_probe_manually():
    # requires conda kilosort env
    from kilosort.io import save_probe

    # for the dual A1 / PPC implant with a single linear probe in A1 (32 channels)
    # and a 4-shank probe in PPC (another 32 channels)

    cn_count = 64  # channels to sort (no ADC/DAC channels)
    channels = np.arange(cn_count)
    kcoords  = np.zeros(cn_count)  # shank numbers
    for i in range(4):
        kcoords[32 + i*8:32 + (i+1)*8] = i+1

    # probes configurations (site coordinates)
    xc, yc = np.zeros(cn_count), np.zeros(cn_count)

    # A1 probe is 1 shank, -5mm ML, Y-spacing is 35um
    xc[0:32] = -5000
    yc[0:32] = np.array([35*i for i in range(32)])

    # PPC probe is 4 shank, 3mm ML, V-shape
    for i in range(4):
        xc[32 + i*8:32 + (i+1)*8] = -4000 + np.array([i*200 + 3*j for j in range(8)])
        yc[32 + i*8:32 + (i+1)*8] = np.abs(np.array([20*j for j in range(8)]) - 20*4)

    probe = {
        'chanMap': channels,
        'xc': xc,
        'yc': yc,
        'kcoords': kcoords,
        'n_chan': cn_count
    }

    save_probe(probe, 'probe.json')


# def oe2kilosort(oe_xml_path):
#     # LEGACY VERSION. FOR SINGLE PROBE IMPLANTS ONLY. NOW IN utils/probe.py
#     # extract probe information from the OpenEphys XML file and 
#     # store in kilosort probe format dictionnary.
#     # Works for Neuropixels OneBox ONLY for now.
#     root = ET.parse(oe_xml_path).getroot()

#     try:
#         onebox = [x for x in root.findall('SIGNALCHAIN')[0].findall('PROCESSOR') if dict(x.items())['name'] == 'OneBox'][0]
#     except IndexError:
#         raise ValueError('OneBox processor not found in the OpenEphys signal chain')

#     # get sampling rate from probe settings - just in case
#     probe_stream = [x for x in onebox.findall('STREAM') if dict(x.items())['name'] == 'ProbeA'][0]
#     sample_rate = float(probe_stream.attrib['sample_rate'])

#     probe = onebox.findall('EDITOR')[0].findall('NP_PROBE')[0]
#     channels = probe.findall('CHANNELS')[0]
#     x_pos = probe.findall('ELECTRODE_XPOS')[0]
#     y_pos = probe.findall('ELECTRODE_YPOS')[0]

#     xc = [float(x) for x in x_pos.attrib.values()]
#     yc = [float(y) for y in y_pos.attrib.values()]
#     kcoords = [int(desc.split(':')[1]) for ch, desc in channels.attrib.items()]
#     ch_count = len(channels.attrib)

#     # get channel map
#     ch_map = [x for x in root.findall('SIGNALCHAIN')[0].findall('PROCESSOR') if dict(x.items())['name'] == 'Channel Map'][0]
#     channels_xml = list(list(ch_map.findall('CUSTOM_PARAMETERS')[0])[0])  # taking the first "stream"

#     # compute mapping indices
#     channel_map_list = [int(ch.attrib['index']) for ch in channels_xml]
#     if 384 in channel_map_list: channel_map_list.remove(384)

#     channels_probe = [int(ch[2:]) for ch in channels.attrib.keys()]
#     idxs_channels = [channels_probe.index(ch) for ch in channel_map_list]

#     return {
#         'chanMap': [x for x in range(ch_count)],
#         'xc': [xc[i] for i in idxs_channels],
#         'yc': [yc[i] for i in idxs_channels],
#         'kcoords': [kcoords[i] for i in idxs_channels],
#         'n_chan': ch_count
#     }

