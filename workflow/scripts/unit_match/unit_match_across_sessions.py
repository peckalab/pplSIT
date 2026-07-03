import json
import os
import re
import sys
import time
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

import UnitMatchPy.assign_unique_id as aid
import UnitMatchPy.bayes_functions as bf
import UnitMatchPy.default_params as default_params
import UnitMatchPy.GUI as gui
import UnitMatchPy.metric_functions as mf
import UnitMatchPy.overlord as ov
import UnitMatchPy.save_utils as su
import UnitMatchPy.utils as util


def log_progress(message):
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{now}] {message}", flush=True)


def log_timed_step(message):
    class TimedStep:
        def __enter__(self):
            self.start = time.perf_counter()
            log_progress(f"START {message}")
            return self

        def __exit__(self, exc_type, exc, tb):
            elapsed = time.perf_counter() - self.start
            status = "FAILED" if exc_type else "DONE"
            log_progress(f"{status} {message} ({elapsed / 60:.1f} min)")
            return False

    return TimedStep()


def normalize_methods(methods):
    if isinstance(methods, str):
        methods = [methods]
    return [str(method).strip().lower() for method in methods]


def selected_method():
    if hasattr(snakemake.wildcards, "method"):
        return str(snakemake.wildcards.method).strip().lower()
    unit_match_config = snakemake.config.get("unit_match", {})
    return normalize_methods(unit_match_config.get("methods", unit_match_config.get("method", "unitmatchpy")))[0]


def add_deepunitmatch_source_path():
    source_path = snakemake.config.get("unit_match", {}).get("deepunitmatch", {}).get("source_path")
    if not source_path:
        return

    source_path = os.path.abspath(os.path.expanduser(str(source_path)))
    if source_path not in sys.path:
        sys.path.insert(0, source_path)


def parse_stream_waveform_path(path):
    path = Path(path)
    if path.parent.name != "RawWaveforms":
        raise ValueError(f"Expected RawWaveforms/*.npy path, got: {path}")
    if path.parents[1].name == "bombcell":
        stream_dir = path.parents[2]
        waveform_source = "bombcell"
    else:
        stream_dir = path.parents[1]
        waveform_source = "unitmatch"
    if stream_dir.parent.name != "kilosort":
        raise ValueError(f"Expected .../kilosort/<stream>/RawWaveforms/*.npy path, got: {path}")
    return {
        "session": stream_dir.parent.parent.name,
        "stream": stream_dir.name,
        "ks_dir": str(stream_dir),
        "waveform_dir": str(path.parent),
        "waveform_path": str(path),
        "waveform_source": waveform_source,
    }


def session_sort_key(session):
    match = re.search(r"\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}", session)
    if match:
        return (0, match.group(0), session)
    return (1, session)


def record_sort_key(record):
    return (*session_sort_key(record["session"]), record["stream"])


def ks_dir_sort_key(ks_dir):
    stream_dir = Path(ks_dir)
    session = stream_dir.parent.parent.name
    return (*session_sort_key(session), stream_dir.name)


def input_waveform_files():
    input_npy_files = [
        str(path)
        for path in snakemake.input.npy_files
        if str(path).endswith(".npy") and Path(str(path)).parent.name == "RawWaveforms"
    ]

    source = snakemake.config.get("unit_match", {}).get("waveform_source", "unitmatch").strip().lower()
    suffix = Path("RawWaveforms") if source == "unitmatch" else Path("bombcell") / "RawWaveforms"
    animal = snakemake.wildcards.animal
    scanned_npy_files = []
    for session in sorted(snakemake.config.get("session_IDs", []), key=session_sort_key):
        if session.split("_")[0] != animal:
            continue
        kilosort_dir = Path(snakemake.config["dst_path"]) / animal / session / "kilosort"
        if not kilosort_dir.is_dir():
            continue
        for stream_dir in sorted(path for path in kilosort_dir.iterdir() if path.is_dir()):
            wave_dir = stream_dir / suffix
            scanned_npy_files.extend(str(path) for path in sorted(wave_dir.glob("Unit*_RawSpikes.npy")))

    npy_files = sorted(set(input_npy_files).union(scanned_npy_files))
    if scanned_npy_files:
        print(
            f"Found {len(scanned_npy_files)} RawWaveforms .npy files by runtime scan "
            f"across configured sessions ({len(input_npy_files)} were listed by Snakemake inputs).",
            flush=True,
        )
    return npy_files


def probe_geometry_signature(probe_json):
    required_keys = ("xc", "yc", "kcoords")
    missing_keys = [key for key in required_keys if key not in probe_json]
    if missing_keys:
        raise ValueError(f"probe.json is missing required geometry key(s): {missing_keys}")

    coords = sorted(
        (float(x), float(y), int(kcoord))
        for x, y, kcoord in zip(probe_json["xc"], probe_json["yc"], probe_json["kcoords"])
    )
    metadata = {
        "probe_serial_number": probe_json.get("probe_serial_number"),
        "headstage_serial_number": probe_json.get("headstage_serial_number"),
        "probe_name": probe_json.get("probe_name"),
        "dock": probe_json.get("dock"),
        "n_chan": probe_json.get("n_chan"),
        "sample_rate": probe_json.get("sample_rate"),
    }
    return json.dumps({"metadata": metadata, "coords": coords}, sort_keys=True)


def exact_probe_signature(probe_json):
    return json.dumps(probe_json, sort_keys=True)


def probe_group_name(group_index):
    return f"probe_geometry_{group_index}"


def warn_if_probe_channel_order_differs(records):
    exact_signatures = defaultdict(list)
    for record in records:
        probe_path = Path(record["ks_dir"]) / "probe.json"
        with open(probe_path, "r") as f:
            probe_json = json.load(f)
        exact_signatures[exact_probe_signature(probe_json)].append(record)

    if len(exact_signatures) <= 1:
        return

    print(
        "Warning: grouping recordings with the same physical probe geometry but "
        f"{len(exact_signatures)} different exact probe.json channel orderings.",
        flush=True,
    )
    for order_index, order_records in enumerate(
        sorted(exact_signatures.values(), key=lambda group: (-len(group), group[0]["session"])),
        start=1,
    ):
        recording_ids = ";".join(f"{record['session']}::{record['stream']}" for record in order_records)
        print(
            f"  channel_order_{order_index}: {len(order_records)} recording(s): {recording_ids}",
            flush=True,
        )


def grouped_recordings_from_inputs(npy_files):
    require_phy_folder = snakemake.config.get("unit_match", {}).get("require_phy_folder", False)
    recording_to_npy_files = {}
    recording_records = {}

    for npy_path in npy_files:
        record = parse_stream_waveform_path(npy_path)
        stream_dir = Path(record["ks_dir"])
        if require_phy_folder and not (stream_dir / ".phy").is_dir():
            continue
        if not (stream_dir / "probe.json").exists():
            print(f"Skipping {stream_dir}: missing probe.json.")
            continue

        recording_key = (record["session"], record["stream"])
        recording_to_npy_files.setdefault(recording_key, []).append(npy_path)
        recording_records[recording_key] = {
            "session": record["session"],
            "stream": record["stream"],
            "ks_dir": record["ks_dir"],
            "waveform_dir": str(Path(record["waveform_path"]).parent),
        }

    groups = {}
    for record in recording_records.values():
        probe_path = Path(record["ks_dir"]) / "probe.json"
        with open(probe_path, "r") as f:
            probe_json = json.load(f)
        probe_key = probe_geometry_signature(probe_json)
        groups.setdefault(probe_key, []).append(record)

    grouped = []
    for _, records in sorted(
        groups.items(),
        key=lambda item: (-len(item[1]), [record_sort_key(r) for r in item[1]]),
    ):
        group_records = sorted(records, key=record_sort_key)
        group_npy_files = []
        for record in group_records:
            group_npy_files.extend(recording_to_npy_files[(record["session"], record["stream"])])
        warn_if_probe_channel_order_differs(group_records)
        if group_records:
            print(
                "Chronological UnitMatch session order: "
                f"{group_records[0]['session']}::{group_records[0]['stream']} -> "
                f"{group_records[-1]['session']}::{group_records[-1]['stream']} "
                f"({len(group_records)} recording(s)).",
                flush=True,
            )
        grouped.append((group_records, group_npy_files))

    return grouped


def load_unitmatch_inputs(group_npy_files):
    records = [parse_stream_waveform_path(path) for path in group_npy_files]
    ks_to_waveform_dir = {}
    for record in records:
        ks_to_waveform_dir.setdefault(record["ks_dir"], record["waveform_dir"])
    ks_dirs = sorted(ks_to_waveform_dir, key=ks_dir_sort_key)
    waveform_dirs = [ks_to_waveform_dir[ks_dir] for ks_dir in ks_dirs]
    log_progress(
        f"Loading UnitMatch inputs from {len(ks_dirs)} recording(s) and "
        f"{len(group_npy_files)} waveform file(s)."
    )
    param = default_params.get_default_param(param={"KS_dirs": ks_dirs})
    wave_paths, unit_label_paths, channel_pos = util.paths_from_KS(
        ks_dirs,
        param=param,
        custom_raw_waveform_paths=waveform_dirs,
    )
    param = util.get_probe_geometry(channel_pos[0], param)
    waveform, session_id, session_switch, within_session, good_units, param = util.load_good_waveforms(
        wave_paths,
        unit_label_paths,
        param,
        good_units_only=snakemake.config["unit_match"]["match_good_units_only"],
    )
    param["good_units"] = good_units
    log_progress(
        f"Loaded waveform array with shape {waveform.shape}; "
        f"{len(np.concatenate(good_units))} unit(s) across {len(ks_dirs)} recording(s)."
    )

    if not np.isfinite(waveform).all():
        n_bad = np.size(waveform) - np.count_nonzero(np.isfinite(waveform))
        print(f"Found {n_bad} non-finite values in waveform; replacing with 0.")
        waveform = np.nan_to_num(waveform, nan=0.0, posinf=0.0, neginf=0.0)

    clus_info = {
        "good_units": good_units,
        "session_switch": session_switch,
        "session_id": session_id,
        "original_ids": np.concatenate(good_units),
    }
    return ks_dirs, waveform, session_id, session_switch, within_session, good_units, channel_pos, param, clus_info


def patch_drift_n_sessions():
    original = mf.drift_n_sessions

    def drift_n_sessions_safe(candidate_pairs, session_switch, avg_centroid, avg_waveform_per_tp, total_score, param):
        try:
            pairs = np.asarray(candidate_pairs)
            if pairs.size == 0 or pairs.ndim < 2:
                n_sessions = len(np.unique(session_switch))
                drifts = np.zeros(n_sessions)
                print("Skipping drift correction: no valid across-session pairs.")
                return drifts, avg_centroid, avg_waveform_per_tp
            return original(
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

    mf.drift_n_sessions = drift_n_sessions_safe


def save_session_mapping(group_records, save_dir):
    session_id_df = pd.DataFrame({
        "recording_id": [f"{record['session']}::{record['stream']}" for record in group_records],
        "session_id": [record["session"] for record in group_records],
        "stream": [record["stream"] for record in group_records],
        "ks_dir": [record["ks_dir"] for record in group_records],
        "session_index": range(len(group_records)),
    })
    session_id_df.to_csv(os.path.join(save_dir, "session_id_mapping.csv"), index=False)
    log_progress(f"Saved session_id_mapping.csv for {len(group_records)} recording(s).")


def save_unitmatch_output(
    save_dir,
    scores_to_include,
    matches,
    output_prob_matrix,
    extracted_wave_properties,
    total_score,
    output_threshold,
    clus_info,
    param,
):
    with log_timed_step(f"assigning unique IDs for {output_prob_matrix.shape[0]} unit(s)"):
        uids = aid.assign_unique_id(output_prob_matrix, param, clus_info)
    with log_timed_step(f"saving UnitMatch output to {save_dir}"):
        su.save_to_output(
            save_dir,
            scores_to_include,
            matches,
            output_prob_matrix,
            extracted_wave_properties["avg_centroid"],
            extracted_wave_properties["avg_waveform"],
            extracted_wave_properties["avg_waveform_per_tp"],
            extracted_wave_properties["max_site"],
            total_score,
            output_threshold,
            clus_info,
            param,
            UIDs=uids,
            matches_curated=None,
            save_match_table=True,
        )


def filter_units_by_index_compat(waveform, session_id, session_switch, good_units, kept_idx, param):
    if hasattr(util, "filter_units_by_index"):
        return util.filter_units_by_index(waveform, session_id, session_switch, good_units, kept_idx, param)

    kept_idx = np.asarray(kept_idx, dtype=int)
    old_session_id = np.asarray(session_id)
    original_ids = np.concatenate(good_units)

    waveform = waveform[kept_idx]
    session_id = old_session_id[kept_idx]
    filtered_ids = original_ids[kept_idx]

    filtered_good_units = []
    session_counts = []
    for session in sorted(np.unique(old_session_id)):
        session_units = filtered_ids[session_id == session]
        filtered_good_units.append(session_units)
        session_counts.append(len(session_units))

    session_switch = np.concatenate(([0], np.cumsum(session_counts)))
    within_session = session_id[:, None] == session_id[None, :]
    param["good_units"] = filtered_good_units
    param["n_units"] = len(filtered_ids)
    return waveform, session_id, session_switch, within_session, filtered_good_units, param


def run_unitmatchpy_for_group(group_npy_files, group_records, save_dir):
    _, waveform, session_id, session_switch, within_session, _, channel_pos, param, clus_info = load_unitmatch_inputs(group_npy_files)

    patch_drift_n_sessions()
    with log_timed_step("UnitMatchPy waveform parameter extraction"):
        extracted_wave_properties = ov.extract_parameters(waveform, channel_pos, clus_info, param)
    with log_timed_step("UnitMatchPy metric scoring"):
        total_score, candidate_pairs, scores_to_include, predictors = ov.extract_metric_scores(
            extracted_wave_properties,
            session_switch,
            within_session,
            param,
            niter=snakemake.config.get("unit_match", {}).get("metrics_niter", 2),
        )

    prior_match = 1 - (param["n_expected_matches"] / param["n_units"] ** 2)
    priors = np.array((prior_match, 1 - prior_match))
    labels = candidate_pairs.astype(int)
    cond = np.unique(labels)
    parameter_kernels = bf.get_parameter_kernels(scores_to_include, labels, cond, param, add_one=1)
    probability = bf.apply_naive_bayes(parameter_kernels, priors, predictors, param, cond)
    output_prob_matrix = probability[:, 1].reshape(param["n_units"], param["n_units"])

    with log_timed_step("UnitMatchPy output evaluation"):
        util.evaluate_output(output_prob_matrix, param, within_session, session_switch, match_threshold=0.75)
    match_threshold = param["match_threshold"]
    output_threshold = np.zeros_like(output_prob_matrix)
    output_threshold[output_prob_matrix > match_threshold] = 1

    if not snakemake.config.get("unit_match", {}).get("skip_gui", False):
        gui.process_info_for_GUI(
            output_prob_matrix,
            match_threshold,
            scores_to_include,
            total_score,
            extracted_wave_properties["amplitude"],
            extracted_wave_properties["spatial_decay"],
            extracted_wave_properties["avg_centroid"],
            extracted_wave_properties["avg_waveform"],
            extracted_wave_properties["avg_waveform_per_tp"],
            extracted_wave_properties["good_wave_idxs"],
            extracted_wave_properties["max_site"],
            extracted_wave_properties["max_site_mean"],
            waveform,
            within_session,
            channel_pos,
            clus_info,
            param,
        )

    save_unitmatch_output(
        save_dir,
        scores_to_include,
        np.argwhere(output_threshold == 1),
        output_prob_matrix,
        extracted_wave_properties,
        total_score,
        output_threshold,
        clus_info,
        param,
    )
    save_session_mapping(group_records, save_dir)


def run_deepunitmatch_for_group(group_npy_files, group_records, save_dir, group_index):
    log_progress(
        f"Preparing DeepUnitMatch for {probe_group_name(group_index)} with "
        f"{len(group_records)} recording(s)."
    )
    add_deepunitmatch_source_path()
    try:
        from DeepUnitMatch.testing import test
        from DeepUnitMatch.utils import helpers, param_fun
    except ImportError as exc:
        raise ImportError(
            "unit_match.methods includes 'deepunitmatch', but DeepUnitMatch is not importable. "
            "Install/rebuild the UnitMatch environment with DeepUnitMatch support, preferably "
            "Python 3.11 or 3.12 plus UnitMatchPy and torch. If DeepUnitMatch is available only "
            "from a source checkout, set unit_match.deepunitmatch.source_path to the checkout's "
            "UnitMatchPy directory."
        ) from exc

    _, waveform, session_id, session_switch, within_session, _, channel_pos, param, clus_info = load_unitmatch_inputs(group_npy_files)
    deep_config = snakemake.config.get("unit_match", {}).get("deepunitmatch", {})
    device = deep_config.get("device", "cuda")
    threshold = float(deep_config.get("threshold", 0.5))
    model_path = deep_config.get("model_path")
    expected_spike_width = deep_config.get("expected_spike_width", 82)
    if expected_spike_width is not None and waveform.shape[1] != int(expected_spike_width):
        raise ValueError(
            "DeepUnitMatch pretrained inference expects raw waveforms with "
            f"spike_width={expected_spike_width}, but loaded waveform has spike_width={waveform.shape[1]}. "
            "Re-run raw waveform extraction with unit_match.spike_width=82 and samples_before=20, "
            "or set unit_match.deepunitmatch.expected_spike_width to match a compatible custom model."
        )
    tmp_root = deep_config.get("tmp_dir")
    if tmp_root is None:
        tmp_root = os.path.join(
            snakemake.config["prj_path"],
            "unit_match",
            snakemake.wildcards.animal,
            "deepunitmatch",
            "tmp",
            probe_group_name(group_index),
        )
    os.makedirs(tmp_root, exist_ok=True)
    log_progress(f"DeepUnitMatch scratch directory: {tmp_root}")

    unit_ids = np.concatenate(param["good_units"]).squeeze()
    with log_timed_step("DeepUnitMatch snippet extraction/preparation"):
        _, _, kept_idx = param_fun.get_snippets(
            waveform,
            channel_pos,
            session_id,
            save_path=tmp_root,
            unit_ids=unit_ids,
            param=param,
        )
    log_progress(f"DeepUnitMatch kept {len(kept_idx)} of {len(waveform)} unit(s) after snippet preparation.")
    if len(kept_idx) < len(waveform):
        log_progress("Filtering UnitMatch inputs to DeepUnitMatch-kept units.")
        waveform, session_id, session_switch, within_session, good_units, param = filter_units_by_index_compat(
            waveform,
            session_id,
            session_switch,
            param["good_units"],
            kept_idx,
            param,
        )
        param["good_units"] = good_units
        clus_info = {
            "good_units": good_units,
            "session_switch": session_switch,
            "session_id": session_id,
            "original_ids": np.concatenate(good_units),
        }

    with log_timed_step(f"loading DeepUnitMatch model on device={device}"):
        if model_path:
            log_progress(f"Using DeepUnitMatch model path: {model_path}")
            model = test.load_trained_model(read_path=model_path, device=device)
        else:
            log_progress("Using DeepUnitMatch default pretrained model.")
            model = test.load_trained_model(device=device)

    data_dir = os.path.join(tmp_root, "processed_waveforms")
    with log_timed_step(f"DeepUnitMatch inference over {len(np.unique(session_id))} session(s)"):
        sim_matrix = test.inference(model, data_dir)
    log_progress(f"DeepUnitMatch similarity matrix shape: {sim_matrix.shape}")

    with log_timed_step("UnitMatch parameter extraction for DeepUnitMatch post-processing"):
        extracted_wave_properties = ov.extract_parameters(waveform, channel_pos, clus_info, param)
    within_session = 1 - (session_id[:, None] == session_id).astype(int)
    sessions = np.unique(session_id)
    probs = np.zeros(sim_matrix.shape)
    distance_matrix = np.zeros(sim_matrix.shape)

    n_pairs = len(sessions) * (len(sessions) - 1) // 2
    pair_count = 0
    pair_start = time.perf_counter()
    log_progress(f"START DeepUnitMatch pairwise Bayesian post-processing for {n_pairs} session pair(s).")
    for r1 in sessions:
        for r2 in sessions:
            if r1 >= r2:
                continue
            pair_count += 1
            if pair_count == 1 or pair_count % 25 == 0 or pair_count == n_pairs:
                elapsed = time.perf_counter() - pair_start
                rate = pair_count / elapsed if elapsed > 0 else 0
                remaining = (n_pairs - pair_count) / rate if rate > 0 else 0
                log_progress(
                    "DeepUnitMatch pairwise post-processing "
                    f"{pair_count}/{n_pairs} session pair(s) "
                    f"(elapsed {elapsed / 60:.1f} min, eta {remaining / 60:.1f} min)."
                )

            mask = np.isin(session_id, [r1, r2])
            sim_mat = sim_matrix[mask][:, mask]
            n_units_r1 = session_switch[r1 + 1] - session_switch[r1]
            n_units_r2 = session_switch[r2 + 1] - session_switch[r2]
            session_switch_pair = np.array([0, n_units_r1, n_units_r1 + n_units_r2])
            indices = np.where(mask)[0]

            df = helpers.create_dataframe(
                [param["good_units"][r1], param["good_units"][r2]],
                sim_mat,
                session_list=[r1, r2],
            )
            matches = test.get_matches(df, sim_mat, session_id[indices], data_dir, dist_thresh=50)

            labels = np.eye(sim_mat.shape[0])
            subsession_id = np.array(
                [r1] * len(param["good_units"][r1]) + [r2] * len(param["good_units"][r2])
            )
            for (rec_ses_1, rec_ses_2), group in matches.groupby(by=["RecSes1", "RecSes2"]):
                as_matrix = group["match"].values.reshape(
                    len(param["good_units"][rec_ses_1]),
                    len(param["good_units"][rec_ses_2]),
                ).astype(int)
                labels[np.ix_(subsession_id == rec_ses_1, subsession_id == rec_ses_2)] = as_matrix

            avg_centroid = extracted_wave_properties["avg_centroid"][:, mask, :]
            avg_waveform_per_tp = extracted_wave_properties["avg_waveform_per_tp"][:, mask, :, :]
            avg_waveform_per_tp = mf.drift_correct_session_pair(
                labels.astype(bool),
                session_switch_pair,
                avg_centroid,
                avg_waveform_per_tp,
                0,
                param,
            )
            avg_waveform_per_tp_flip = mf.flip_dim(avg_waveform_per_tp, param, np.sum(mask))
            euclid_dist = mf.get_Euclidean_dist(avg_waveform_per_tp_flip, param, np.sum(mask))
            centroid_dist, _ = mf.centroid_metrics(euclid_dist, param)
            scores_to_include = {
                "similarity": sim_mat,
                "distance": centroid_dist,
            }

            n_units = int(np.sqrt(len(df)))
            priors = np.array([1 - 2 / n_units, 2 / n_units])
            cond = np.unique(labels)
            parameter_kernels = bf.get_parameter_kernels(scores_to_include, labels, cond, param)
            predictors = np.stack([scores for scores in scores_to_include.values()], axis=2)
            probability = bf.apply_naive_bayes(parameter_kernels, priors, predictors, param, cond)
            probs[np.ix_(mask, mask)] = probability[:, 1].reshape(n_units, n_units)
            distance_matrix[np.ix_(mask, mask)] = centroid_dist
    log_progress(f"DONE DeepUnitMatch pairwise Bayesian post-processing ({(time.perf_counter() - pair_start) / 60:.1f} min)")

    with log_timed_step("DeepUnitMatch output evaluation"):
        util.evaluate_output(probs, param, within_session, session_switch, match_threshold=threshold)
    with log_timed_step("DeepUnitMatch directional filtering"):
        final_matches = test.directional_filter(probs, session_id, threshold)
    log_progress(f"Found {np.sum(final_matches) // 2} DeepUnitMatch matches at threshold {threshold}.")

    save_unitmatch_output(
        save_dir,
        {"distance": distance_matrix},
        np.argwhere(final_matches),
        probs,
        extracted_wave_properties,
        distance_matrix,
        final_matches,
        clus_info,
        param,
    )
    save_session_mapping(group_records, save_dir)


def run_method_for_group(method, group_npy_files, group_records, save_dir, group_index):
    if method == "unitmatchpy":
        run_unitmatchpy_for_group(group_npy_files, group_records, save_dir)
    elif method == "deepunitmatch":
        run_deepunitmatch_for_group(group_npy_files, group_records, save_dir, group_index)
    else:
        raise ValueError(f"Unsupported unit_match method: {method!r}")


method = selected_method()
configured_methods = normalize_methods(
    snakemake.config.get("unit_match", {}).get(
        "methods",
        snakemake.config.get("unit_match", {}).get("method", "unitmatchpy"),
    )
)
if method not in configured_methods:
    raise ValueError(f"Method wildcard {method!r} is not listed in unit_match.methods={configured_methods!r}.")

npy_files = input_waveform_files()
grouped_recordings = grouped_recordings_from_inputs(npy_files)

group_results = []
for group_index, (group_records, group_npy_files) in enumerate(grouped_recordings, start=1):
    if not group_npy_files:
        recording_ids = [f"{record['session']}::{record['stream']}" for record in group_records]
        print(f"Skipping {probe_group_name(group_index)}: no RawWaveforms inputs for recordings {recording_ids}.")
        continue

    save_dir = os.path.join(
        snakemake.config["prj_path"],
        "unit_match",
        snakemake.wildcards.animal,
        method,
        probe_group_name(group_index),
    )
    os.makedirs(save_dir, exist_ok=True)

    run_method_for_group(method, group_npy_files, group_records, save_dir, group_index)

    group_results.append({
        "method": method,
        "group": group_index,
        "n_recordings": len(group_records),
        "recordings": ";".join(f"{record['session']}::{record['stream']}" for record in group_records),
        "sessions": ";".join(sorted({record["session"] for record in group_records})),
        "streams": ";".join(sorted({record["stream"] for record in group_records})),
        "save_dir": save_dir,
    })

summary_path = snakemake.output.group_summary
os.makedirs(os.path.dirname(summary_path), exist_ok=True)
pd.DataFrame(group_results).to_csv(summary_path, index=False)
