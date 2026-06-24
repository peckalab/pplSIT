import UnitMatchPy.bayes_functions as bf
import UnitMatchPy.utils as util
import UnitMatchPy.overlord as ov
import numpy as np
import matplotlib.pyplot as plt
import UnitMatchPy.save_utils as su
import UnitMatchPy.GUI as gui
import UnitMatchPy.assign_unique_id as aid
import UnitMatchPy.default_params as default_params
import UnitMatchPy.metric_functions as mf
from pathlib import Path
import os
import pandas as pd
import json


def parse_stream_waveform_path(path):
    path = Path(path)
    if path.parent.name != "RawWaveforms":
        raise ValueError(f"Expected RawWaveforms/*.npy path, got: {path}")
    stream_dir = path.parents[1]
    if stream_dir.parent.name != "kilosort":
        raise ValueError(f"Expected .../kilosort/<stream>/RawWaveforms/*.npy path, got: {path}")
    return {
        "session": stream_dir.parent.parent.name,
        "stream": stream_dir.name,
        "ks_dir": str(stream_dir),
        "waveform_path": str(path),
    }


def group_recordings_by_probe(animal_id, base_path, require_phy_folder=False):
    animal_path = Path(base_path) / animal_id
    if not animal_path.exists():
        raise FileNotFoundError(f"Animal path not found: {animal_path}")

    groups = {}
    for session_dir in sorted(p for p in animal_path.iterdir() if p.is_dir()):
        kilosort_dir = session_dir / "kilosort"
        if not kilosort_dir.is_dir():
            continue

        for stream_dir in sorted(p for p in kilosort_dir.iterdir() if p.is_dir()):
            if require_phy_folder and not (stream_dir / ".phy").is_dir():
                continue

            probe_path = stream_dir / "probe.json"
            if not probe_path.exists():
                continue

            with open(probe_path, "r") as f:
                content = json.load(f)

            key = json.dumps(content, sort_keys=True)
            groups.setdefault(key, []).append({
                "session": session_dir.name,
                "stream": stream_dir.name,
                "ks_dir": str(stream_dir),
            })

    return sorted(
        groups.items(),
        key=lambda item: (-len(item[1]), [(r["session"], r["stream"]) for r in item[1]])
    )


def run_unit_match_for_group(group_npy_files, group_records, save_dir):
    # Get default parameters, can add your own before or after!
    param = default_params.get_default_param()

    # also get KS_dirs from the raw_waveforms paths
    KS_dirs = [parse_stream_waveform_path(p)["ks_dir"] for p in group_npy_files]
    KS_dirs = sorted(list(set(KS_dirs)))

    param['KS_dirs'] = KS_dirs
    wave_paths, unit_label_paths, channel_pos = util.paths_from_KS(KS_dirs)
    param = util.get_probe_geometry(channel_pos[0], param)

    # STEP 0 -- data preparation
    # Read in data and select the good units and exact metadata
    # if snakemake.config['unit_match']['match_good_units_only'] is False, all units will be "good" units
    waveform, session_id, session_switch, within_session, good_units, param = util.load_good_waveforms(
        wave_paths,
        unit_label_paths,
        param,
        good_units_only=snakemake.config['unit_match']['match_good_units_only'],
    )

    # Ensure waveform is finite to avoid scipy.signal.detrend failures
    if not np.isfinite(waveform).all():
        n_bad = np.size(waveform) - np.count_nonzero(np.isfinite(waveform))
        print(f"Found {n_bad} non-finite values in waveform; replacing with 0.")
        waveform = np.nan_to_num(waveform, nan=0.0, posinf=0.0, neginf=0.0)

    # param['peak_loc'] = #may need to set as a value if the peak location is NOT ~ half the spike width

    # Create clus_info, contains all unit id/session related info
    clus_info = {
        'good_units': good_units,
        'session_switch': session_switch,
        'session_id': session_id,
        'original_ids': np.concatenate(good_units),
    }

    # STEP 1
    # Extract parameters from waveform
    extracted_wave_properties = ov.extract_parameters(waveform, channel_pos, clus_info, param)

    # Guard against empty/1D candidate pairs causing drift correction failure
    _orig_drift_n_sessions = mf.drift_n_sessions

    def _drift_n_sessions_safe(candidate_pairs, session_switch, avg_centroid, avg_waveform_per_tp, total_score, param):
        try:
            pairs = np.asarray(candidate_pairs)
            if pairs.size == 0 or pairs.ndim < 2:
                n_sessions = len(np.unique(session_switch))
                drifts = np.zeros(n_sessions)
                print("Skipping drift correction: no valid across-session pairs.")
                return drifts, avg_centroid, avg_waveform_per_tp
            return _orig_drift_n_sessions(
                candidate_pairs,
                session_switch,
                avg_centroid,
                avg_waveform_per_tp,
                total_score,
                param,
            )
        except Exception as exc:
            n_sessions = len(np.unique(session_switch))
            drifts = np.zeros(n_sessions)
            print(f"Skipping drift correction due to error: {exc}")
            return drifts, avg_centroid, avg_waveform_per_tp

    mf.drift_n_sessions = _drift_n_sessions_safe

    # STEP 2, 3, 4
    # Extract metric scores
    total_score, candidate_pairs, scores_to_include, predictors = ov.extract_metric_scores(
        extracted_wave_properties,
        session_switch,
        within_session,
        param,
        niter=2,
    )

    # STEP 5
    # Probability analysis
    # Get prior probability of being a match
    prior_match = 1 - (param['n_expected_matches'] / param['n_units']**2)  # freedom of choose in prior prob
    priors = np.array((prior_match, 1 - prior_match))

    # Construct distributions (kernels) for Naive Bayes Classifier
    labels = candidate_pairs.astype(int)
    cond = np.unique(labels)
    score_vector = param['score_vector']
    parameter_kernels = np.full((len(score_vector), len(scores_to_include), len(cond)), np.nan)

    parameter_kernels = bf.get_parameter_kernels(scores_to_include, labels, cond, param, add_one=1)

    # Get probability of each pair of being a match
    probability = bf.apply_naive_bayes(parameter_kernels, priors, predictors, param, cond)

    output_prob_matrix = probability[:, 1].reshape(param['n_units'], param['n_units'])

    util.evaluate_output(output_prob_matrix, param, within_session, session_switch, match_threshold=0.75)

    match_threshold = param['match_threshold']
    # match_threshold = try different values here!

    output_threshold = np.zeros_like(output_prob_matrix)
    output_threshold[output_prob_matrix > match_threshold] = 1

    # STEP 6
    # Format data for GUI
    amplitude = extracted_wave_properties['amplitude']
    spatial_decay = extracted_wave_properties['spatial_decay']
    avg_centroid = extracted_wave_properties['avg_centroid']
    avg_waveform = extracted_wave_properties['avg_waveform']
    avg_waveform_per_tp = extracted_wave_properties['avg_waveform_per_tp']
    wave_idx = extracted_wave_properties['good_wave_idxs']
    max_site = extracted_wave_properties['max_site']
    max_site_mean = extracted_wave_properties['max_site_mean']
    gui.process_info_for_GUI(
        output_prob_matrix,
        match_threshold,
        scores_to_include,
        total_score,
        amplitude,
        spatial_decay,
        avg_centroid,
        avg_waveform,
        avg_waveform_per_tp,
        wave_idx,
        max_site,
        max_site_mean,
        waveform,
        within_session,
        channel_pos,
        clus_info,
        param,
    )

    matches = np.argwhere(output_threshold == 1)
    UIDs = aid.assign_unique_id(output_prob_matrix, param, clus_info)

    # NOTE - change to matches to matches_curated if done manual curation with the GUI
    su.save_to_output(
        save_dir,
        scores_to_include,
        matches,  # matches_curated
        output_prob_matrix,
        avg_centroid,
        avg_waveform,
        avg_waveform_per_tp,
        max_site,
        total_score,
        output_threshold,
        clus_info,
        param,
        UIDs=UIDs,
        matches_curated=None,
        save_match_table=True,
    )

    # save also the correspondence between session ids and their corresponding index
    recording_ids = [f"{record['session']}::{record['stream']}" for record in group_records]
    session_id_df = pd.DataFrame({
        'recording_id': recording_ids,
        'session_id': [record['session'] for record in group_records],
        'stream': [record['stream'] for record in group_records],
        'ks_dir': [record['ks_dir'] for record in group_records],
        'session_index': range(len(group_records)),
    })
    session_id_df.to_csv(os.path.join(save_dir, 'session_id_mapping.csv'), index=False)


method = snakemake.config.get('unit_match', {}).get('method', 'unitmatchpy')
if method != 'unitmatchpy':
    raise NotImplementedError(
        f"unit_match.method={method!r} is not implemented yet. "
        "The stream-aware adapter boundary is in place; use 'unitmatchpy' until DeepUnitMatch is wired."
    )

# snakemake.input.npy_files contains RawWaveforms/*.npy plus label/spike rerun triggers.
npy_files = [
    str(path) for path in snakemake.input.npy_files
    if str(path).endswith('.npy') and Path(str(path)).parent.name == 'RawWaveforms'
]
animal_id = snakemake.wildcards.animal
base_path = snakemake.config.get('dst_path', '/mnt/nevermind.data-share/ag-grothe/AG_Pecka/data/processed')
require_phy_folder = snakemake.config.get('unit_match', {}).get('require_phy_folder', False)

grouped_recordings = group_recordings_by_probe(animal_id, base_path, require_phy_folder=require_phy_folder)

recording_to_npy_files = {}
for npy_path in npy_files:
    record = parse_stream_waveform_path(npy_path)
    recording_key = (record["session"], record["stream"])
    recording_to_npy_files.setdefault(recording_key, []).append(npy_path)

group_results = []
for group_index, (_, records) in enumerate(grouped_recordings, start=1):
    group_records = sorted(records, key=lambda record: (record["session"], record["stream"]))
    group_npy_files = []
    for record in group_records:
        group_npy_files.extend(recording_to_npy_files.get((record["session"], record["stream"]), []))

    if not group_npy_files:
        recording_ids = [f"{record['session']}::{record['stream']}" for record in group_records]
        print(f"Skipping probe_json_{group_index}: no RawWaveforms inputs for recordings {recording_ids}.")
        continue

    save_dir = os.path.join(snakemake.config['prj_path'], 'unit_match', animal_id, f'probe_json_{group_index}')
    os.makedirs(save_dir, exist_ok=True)

    run_unit_match_for_group(group_npy_files, group_records, save_dir)

    group_results.append({
        'group': group_index,
        'n_recordings': len(group_records),
        'recordings': ';'.join(f"{record['session']}::{record['stream']}" for record in group_records),
        'sessions': ';'.join(sorted(set(record['session'] for record in group_records))),
        'streams': ';'.join(sorted(set(record['stream'] for record in group_records))),
        'save_dir': save_dir,
    })

# Write a summary at the expected snakemake output path
summary_path = snakemake.output.group_summary
summary_df = pd.DataFrame(group_results)
summary_df.to_csv(summary_path, index=False)
