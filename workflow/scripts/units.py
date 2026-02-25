import os
import sys
import json
import h5py
import numpy as np
import scipy.ndimage as ndi

# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.neurosuite import load_clu_res, XMLHero
from utils.kilosort import load_ks_units_before, load_ks_units_after
from utils.spiketrain import instantaneous_rate, spike_idxs
from utils.hdf import create_dataset, H5NAMES
from utils.spatial import place_field_2D, map_stats, get_field_patches
from utils.spatial import bins2meters, cart2pol


# unit metrics to compute
metric_names = (
    H5NAMES.o_maps, H5NAMES.f_maps, H5NAMES.sparsity, H5NAMES.selectivity,
    H5NAMES.spat_info, H5NAMES.peak_FR, H5NAMES.spat_info_ns, H5NAMES.peak_FR_ns,
    H5NAMES.f_patches, H5NAMES.f_COM,
    H5NAMES.pfr_center, H5NAMES.occ_info, H5NAMES.o_patches, H5NAMES.o_COM
)


def _flatten_snakemake_inputs(inp):
    """
    Robustly flatten snakemake.input into a list of file paths (strings).
    Handles InputFiles, named inputs, lists/tuples, and dict-like inputs.
    """
    out = []

    # If it's dict-like (has keys), flatten its values
    if hasattr(inp, "keys") and hasattr(inp, "__getitem__") and not isinstance(inp, (str, bytes)):
        try:
            for k in inp.keys():
                v = inp[k]
                if v is None:
                    continue
                if isinstance(v, (list, tuple)):
                    out.extend([str(x) for x in v if x is not None])
                else:
                    out.append(str(v))
            return out
        except Exception:
            # fall back to iterable flatten below
            pass

    # Otherwise treat as iterable (InputFiles behaves like this)
    for v in inp:
        if v is None:
            continue
        if isinstance(v, (list, tuple)):
            out.extend([str(x) for x in v if x is not None])
        else:
            out.append(str(v))

    return out


def _find_kilosort_marker(inputs):
    """
    Find .../<session>/kilosort/.STAGED in the input list.
    """
    for p in inputs:
        if os.path.basename(p) == ".STAGED" and os.path.basename(os.path.dirname(p)) == "kilosort":
            return p
    # fallback: any .STAGED under kilosort
    for p in inputs:
        if os.path.basename(p) == ".STAGED" and "kilosort" in p.split(os.sep):
            return p
    return None


def _read_streams(marker_path):
    with open(marker_path, "r") as f:
        return [ln.strip() for ln in f if ln.strip()]


def _electrode_offset_for_stream(stream_idx, n_shanks_per_stream=4):
    return stream_idx * n_shanks_per_stream


def _load_timeline(meta_h5_path):
    with h5py.File(meta_h5_path, "r") as f:
        tl = np.array(f["processed"]["timeline"])  # time, X, Y, speed, HD, trials, sounds
    return tl


def _compute_and_store_unit_metrics(out_h5, unit_name, tl, s_times, spiketrain_raw, sampling_rate,
                                   positions=None, unit_info=None, electrode_idx=None, unit_idx=None):
    """
    Writes spike data + spatial metrics for one unit into units.h5 using create_dataset.
    """

    # instant rate + spike indices on timeline
    i_rate = instantaneous_rate(s_times, tl[:, 0])
    s_idxs = spike_idxs(s_times, tl[:, 0])

    create_dataset(out_h5, unit_name, H5NAMES.spike_times, s_times)
    create_dataset(out_h5, unit_name, H5NAMES.inst_rate, i_rate)
    create_dataset(out_h5, unit_name, H5NAMES.spike_idxs, s_idxs)

    # mean firing metrics (keep same semantics as old script)
    mean_rate = len(spiketrain_raw) / (tl[-1][0] - tl[0][0])
    create_dataset(out_h5, unit_name, H5NAMES.mfr, mean_rate)

    # NOTE: old code used spiketrain_raw (samples idxs for KS, sample counts for NS) for robust_rate.
    # We keep that behavior to avoid changing downstream expectations, even if it's imperfect for KS.
    robust_rate = 1 / np.median(np.diff(spiketrain_raw)) if len(spiketrain_raw) > 2 else 0.0
    create_dataset(out_h5, unit_name, H5NAMES.mfr_robust, robust_rate)

    if positions is not None and electrode_idx is not None and unit_idx is not None:
        create_dataset(out_h5, unit_name, H5NAMES.anat_pos, positions[electrode_idx][unit_idx])

    if unit_info is not None and electrode_idx is not None and unit_idx is not None:
        create_dataset(out_h5, unit_name, H5NAMES.kilosort, unit_info[electrode_idx][unit_idx])

    # spatial metrics
    xy_range = [-0.5, 0.5, -0.5, 0.5]  # fixed for cross-comparisons
    bin_size = 0.02  # 2 cm
    s_rate_pos = round(1.0 / np.diff(tl[:, 0]).mean())

    # whole session
    unit_pos = tl[s_idxs][:, 1:3]
    traj_pos = tl[:, 1:3]
    o_map, s1_map, s2_map, f_map = place_field_2D(traj_pos, unit_pos, s_rate_pos, bin_size=bin_size, xy_range=xy_range)
    sparsity, selectivity, spat_info, peak_FR = map_stats(f_map, o_map)

    # no stimulus only
    idxs_tl_nostim = np.where(tl[:, 6] == 0)[0]
    unit_pos_ns = tl[np.intersect1d(s_idxs, idxs_tl_nostim)][:, 1:3]
    traj_pos_ns = tl[idxs_tl_nostim][:, 1:3]
    o_map_ns, s1_map_ns, s2_map_ns, f_map_ns = place_field_2D(traj_pos_ns, unit_pos_ns, s_rate_pos, bin_size=bin_size, xy_range=xy_range)
    _, _, spat_info_ns, peak_FR_ns = map_stats(f_map_ns, o_map_ns)

    # place field metrics
    patches = get_field_patches(f_map)
    if f_map.max() == 0:
        f_COM_rho, f_COM_phi, pfr_rho, pfr_phi = 0, 0, 0, 0
    else:
        x_in_b, y_in_b = ndi.center_of_mass(f_map, labels=patches, index=1)
        f_COM_rho, f_COM_phi = cart2pol(*bins2meters(x_in_b, y_in_b, xy_range))
        x, y = np.where(f_map == np.max(f_map))
        pfr_rho, pfr_phi = cart2pol(*bins2meters(x[0], y[0], xy_range))

    # occupancy metrics
    _, _, occ_info, _ = map_stats(o_map, o_map)
    o_patches = get_field_patches(o_map)
    x, y = ndi.center_of_mass(o_map, labels=o_patches, index=1)
    o_COM_rho, o_COM_phi = cart2pol(*bins2meters(x, y, xy_range))

    # store metrics in the same order as metric_names
    metrics = [
        o_map, f_map, sparsity, selectivity, spat_info, peak_FR, spat_info_ns, peak_FR_ns,
        patches, (f_COM_rho, f_COM_phi), (pfr_rho, pfr_phi), occ_info, o_patches, (o_COM_rho, o_COM_phi)
    ]
    for i, ds in enumerate(metrics):
        create_dataset(out_h5, unit_name, metric_names[i], np.array(ds))


def _process_neurosuite(sorted_data_path, tl, out_h5_path):
    # neurosuite: read from XML
    xml_files = [f for f in os.listdir(sorted_data_path) if f.find(".xml") > 0]
    if len(xml_files) == 0:
        raise ValueError("Need XML settings file to store spiketrains sorted by Neurosuite")

    neurosuite_settings_file = os.path.join(sorted_data_path, xml_files[0])
    sampling_rate = XMLHero(neurosuite_settings_file).get_sampling_rate()
    timestamps_path = None  # neurosuite doesn't use timestamps.npy

    # spikes are in samples, not seconds
    units = load_clu_res(sorted_data_path)

    for electrode_idx in units.keys():
        unit_idxs = units[electrode_idx]
        for unit_idx, spiketrain in unit_idxs.items():
            unit_name = f"{electrode_idx}-{unit_idx}"

            # convert to seconds using sampling rate
            s_times = spiketrain / sampling_rate

            _compute_and_store_unit_metrics(
                out_h5=out_h5_path,
                unit_name=unit_name,
                tl=tl,
                s_times=s_times,
                spiketrain_raw=spiketrain,
                sampling_rate=sampling_rate,
                positions=None,
                unit_info=None,
                electrode_idx=None,
                unit_idx=None
            )


def _process_kilosort_stream(stream_name, stream_folder, tl, out_h5_path, animal, session, config, electrode_offset=0):
    # settings + probe
    kilosort_settings_file = os.path.join(stream_folder, "settings.json")
    clu_info_file = os.path.join(stream_folder, "cluster_info.tsv")

    if not os.path.exists(kilosort_settings_file):
        raise FileNotFoundError(f"Missing settings.json for stream {stream_name}: {kilosort_settings_file}")

    with open(kilosort_settings_file, "r") as json_file:
        sampling_rate = json.load(json_file)["fs"]

    # timestamps: prefer staged per-stream timestamps.npy
    timestamps_path = os.path.join(stream_folder, "timestamps.npy")
    if not os.path.exists(timestamps_path):
        # fallback to raw layout (old behavior)
        timestamps_path = os.path.join(config["src_path"], animal, session, "timestamps.npy")
        if not os.path.exists(timestamps_path):
            raise FileNotFoundError(
                f"No timestamps.npy found for stream {stream_name}. "
                f"Looked in {os.path.join(stream_folder, 'timestamps.npy')} and {timestamps_path}"
            )

    # load units
    positions, unit_info = None, None
    if os.path.exists(clu_info_file):
        units, unit_info = load_ks_units_after(stream_folder)
    else:
        units, positions = load_ks_units_before(stream_folder)

    # write units for this stream
    timestamps = np.load(timestamps_path)

    for electrode_idx in units.keys():
        unit_idxs = units[electrode_idx]
        for unit_idx, spiketrain in unit_idxs.items():
            # prefix to avoid collisions across streams
            #unit_name = f"{stream_name}:{electrode_idx}-{unit_idx}"
            global_electrode_idx = int(electrode_idx) + int(electrode_offset)
            unit_name = f"{global_electrode_idx}-{unit_idx}"

            # kilosort: spiketrain are indices into timestamps (per your loader)
            s_times = timestamps[spiketrain] - timestamps[0]

            _compute_and_store_unit_metrics(
                out_h5=out_h5_path,
                unit_name=unit_name,
                tl=tl,
                s_times=s_times,
                spiketrain_raw=spiketrain,
                sampling_rate=sampling_rate,
                positions=positions,
                unit_info=unit_info,
                electrode_idx=electrode_idx,
                unit_idx=unit_idx
            )


# load timeline once
tl = _load_timeline(snakemake.input[0])

source = snakemake.config["units"]["source"]

if source == "neurosuite":
    # Expect second input to be something like neurosuite.ready inside the sorted folder
    inputs = _flatten_snakemake_inputs(snakemake.input)
    if len(inputs) < 2:
        raise ValueError("Neurosuite mode expects at least 2 inputs: meta.h5 and neurosuite.ready (or similar).")
    sorted_data_path = os.path.dirname(inputs[1])
    _process_neurosuite(sorted_data_path, tl, snakemake.output[0])

else:
    # kilosort multi-stream
    inputs = _flatten_snakemake_inputs(snakemake.input)
    marker = _find_kilosort_marker(inputs)
    if marker is None:
        raise ValueError(
            "Kilosort mode expects kilosort/.STAGED marker in inputs "
            "(e.g. processed/<animal>/<session>/kilosort/.STAGED)."
        )

    streams = _read_streams(marker)
    ks_root = os.path.dirname(marker)  # .../<session>/kilosort

    animal = snakemake.params["animal"]
    session = snakemake.params["session"]

    for stream_idx, stream in enumerate(streams):
        stream_folder = os.path.join(ks_root, stream)

        ready_path = os.path.join(stream_folder, "kilosort.ready")
        if not os.path.exists(ready_path):
            raise FileNotFoundError(f"Missing per-stream kilosort.ready: {ready_path}")

        _process_kilosort_stream(
            stream_name=stream,
            stream_folder=stream_folder,
            tl=tl,
            out_h5_path=snakemake.output[0],
            animal=animal,
            session=session,
            config=snakemake.config,
            electrode_offset=_electrode_offset_for_stream(stream_idx, n_shanks_per_stream=4),
        )



# import os, sys
# import h5py
# import json
# import numpy as np
# import scipy.ndimage as ndi

# # import util functions from utils module
# parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
# sys.path.append(os.getcwd())
# sys.path.append(parent_dir)

# from utils.neurosuite import load_clu_res, XMLHero
# from utils.kilosort import load_ks_units_before, load_ks_units_after
# from utils.spiketrain import instantaneous_rate, spike_idxs
# from utils.hdf import create_dataset, H5NAMES
# from utils.spatial import place_field_2D, map_stats, get_field_patches
# from utils.spatial import bins2meters, cart2pol


# # unit metrics to compute
# metric_names = (H5NAMES.o_maps, H5NAMES.f_maps, H5NAMES.sparsity, H5NAMES.selectivity, \
#                 H5NAMES.spat_info, H5NAMES.peak_FR, H5NAMES.spat_info_ns, H5NAMES.peak_FR_ns, \
#                 H5NAMES.f_patches, H5NAMES.f_COM, \
#                 H5NAMES.pfr_center, H5NAMES.occ_info, H5NAMES.o_patches, H5NAMES.o_COM)


# # loading spike data
# sorted_data_path = os.path.dirname(snakemake.input[1])

# positions, unit_info = None, None
# if snakemake.config['units']['source'] == 'neurosuite':
#     # neurosuite: read from XML
#     xml_files = [f for f in os.listdir(sorted_data_path) if f.find('.xml') > 0]
#     if len(xml_files) == 0:
#         raise ValueError('Need XML settings file to store spiketrains sorted by Neurosuite')

#     neurosuite_settings_file = os.path.join(sorted_data_path, xml_files[0])
#     sampling_rate = XMLHero(neurosuite_settings_file).get_sampling_rate()
#     timestamps_path = None

#     # loading unit data from .clu / .res
#     units = load_clu_res(sorted_data_path)  # spikes are in samples, not seconds

# else:
#     # kilosort: read from settings.json
#     kilosort_settings_file = os.path.join(sorted_data_path, 'settings.json')
#     probe_file = os.path.join(sorted_data_path, 'probe.json')
#     clu_info_file = os.path.join(sorted_data_path, 'cluster_info.tsv')
#     timestamps_path = os.path.join(snakemake.config['src_path'], snakemake.params['animal'], snakemake.params['session'], 'timestamps.npy')
#     with open(kilosort_settings_file, 'r') as json_file:
#         sampling_rate = json.load(json_file)['fs']
#     with open(probe_file, 'r') as json_file:
#         probe = json.load(json_file)

#     # loading unit data from kilosort
#     if os.path.exists(clu_info_file):
#         units, unit_info = load_ks_units_after(sorted_data_path)
#     else:
#         units, positions = load_ks_units_before(sorted_data_path)


# # loading timeline
# with h5py.File(snakemake.input[0], 'r') as f:
#     tl = np.array(f['processed']['timeline'])  # time, X, Y, speed, HD, trials, sounds
#     run_idxs = np.where(tl[:, 3] > 0.04)[0]


# # writing spike in our formats with metrics
# for electrode_idx in units.keys():
#     unit_idxs = units[electrode_idx]

#     for unit_idx, spiketrain in unit_idxs.items():
#         unit_name = '%s-%s' % (electrode_idx, unit_idx)

#         # 3 main ways to store a spiketrain
#         if timestamps_path is not None:
#             # kilosort: read from timestamps
#             timestamps = np.load(timestamps_path)
#             s_times = timestamps[spiketrain] - timestamps[0] # timestamps in seconds
#         else:
#             # neurosuite: assume sampling rate is known
#             s_times = spiketrain / sampling_rate  # spike times in seconds
        
#         i_rate  = instantaneous_rate(s_times, tl[:, 0])  # instantaneous rate, sampling to timeline
#         s_idxs  = spike_idxs(s_times, tl[:, 0])  # timeline indices of individual spikes

#         create_dataset(snakemake.output[0], unit_name, H5NAMES.spike_times, s_times)
#         create_dataset(snakemake.output[0], unit_name, H5NAMES.inst_rate, i_rate)
#         create_dataset(snakemake.output[0], unit_name, H5NAMES.spike_idxs, s_idxs)

#         # mean firing
#         mean_rate = len(spiketrain) / (tl[-1][0] - tl[0][0])  # mean firing rate - spike count / time
#         create_dataset(snakemake.output[0], unit_name, H5NAMES.mfr, mean_rate)

#         robust_rate = 1 / np.median(np.diff(spiketrain))  #  mean firing rate - 1 / median(ISI)
#         create_dataset(snakemake.output[0], unit_name, H5NAMES.mfr_robust, robust_rate)

#         if positions is not None:
#             create_dataset(snakemake.output[0], unit_name, H5NAMES.anat_pos, positions[electrode_idx][unit_idx])
#         if unit_info is not None:
#             create_dataset(snakemake.output[0], unit_name, H5NAMES.kilosort, unit_info[electrode_idx][unit_idx])

#         # spatial metrics
#         xy_range = [-0.5, 0.5, -0.5, 0.5]  # make fixed for cross-comparisons
#         bin_size = 0.02  # 2 cm
#         s_rate_pos = round(1.0 / np.diff(tl[:, 0]).mean())

#         # keep only spiking when running? > 4cm/s
#         #s_idxs = np.intersect1d(s_idxs, run_idxs)

#         # compute 2D maps: occupancy and firing rate (place fields) for 
#         # a) the whole session
#         unit_pos = tl[s_idxs][:, 1:3]
#         traj_pos = tl[:, 1:3]
#         #xy_range = [tl[:, 1].min(), tl[:, 1].max(), tl[:, 2].min(), tl[:, 2].max()]
#         o_map, s1_map, s2_map, f_map = place_field_2D(traj_pos, unit_pos, s_rate_pos, bin_size=bin_size, xy_range=xy_range)
#         sparsity, selectivity, spat_info, peak_FR = map_stats(f_map, o_map)

#         # b) no stimulus periods only
#         idxs_tl_nostim = np.where(tl[:, 6] == 0)[0]
#         unit_pos = tl[np.intersect1d(s_idxs, idxs_tl_nostim)][:, 1:3]
#         traj_pos = tl[idxs_tl_nostim][:, 1:3]
#         #xy_range = [tl[:, 1].min(), tl[:, 1].max(), tl[:, 2].min(), tl[:, 2].max()]
#         o_map_ns, s1_map_ns, s2_map_ns, f_map_ns = place_field_2D(traj_pos, unit_pos, s_rate_pos, bin_size=bin_size, xy_range=xy_range)
#         _, _, spat_info_ns, peak_FR_ns = map_stats(f_map_ns, o_map_ns)

#         # place field metrics
#         patches = get_field_patches(f_map)  # 2D matrix, patches labeled according to the size
#         #f_sizes = np.bincount(patches.flat)[1:]  # 1D array of field sizes, sorted
#         if f_map.max() == 0:
#             f_COM_rho, f_COM_phi, pfr_rho, pfr_phi = 0, 0, 0, 0
#         else:
#             x_in_b, y_in_b = ndi.center_of_mass(f_map, labels=patches, index=1)  # largest field COM, in bins
#             f_COM_rho, f_COM_phi = cart2pol(*bins2meters(x_in_b, y_in_b, xy_range))  # largest field COM, in polar coords.
#             x, y = np.where(f_map == np.max(f_map))  # location of the peak unit firing, in bins
#             pfr_rho, pfr_phi = cart2pol(*bins2meters(x[0], y[0], xy_range))  # location of the peak unit firing, in polar
        
#         # same for occupancy
#         _, _, occ_info, _ = map_stats(o_map, o_map)
#         o_patches = get_field_patches(o_map)  # 2D matrix, patches labeled according to the size
#         x, y = ndi.center_of_mass(o_map, labels=o_patches, index=1)  # largest field COM, in bins
#         o_COM_rho, o_COM_phi = cart2pol(*bins2meters(x, y, xy_range))     # largest field COM, in polar coords.

#         # iterate over metrics, order should match metric_names defined above
#         for i, ds in enumerate([o_map, f_map, sparsity, selectivity, spat_info, peak_FR, spat_info_ns, peak_FR_ns, \
#             patches, (f_COM_rho, f_COM_phi), (pfr_rho, pfr_phi), occ_info, \
#             o_patches, (o_COM_rho, o_COM_phi)]):
#             create_dataset(snakemake.output[0], unit_name, metric_names[i], np.array(ds))