import UnitMatchPy.bayes_functions as bf
import UnitMatchPy.utils as util
import UnitMatchPy.overlord as ov
import numpy as np
import matplotlib.pyplot as plt
import UnitMatchPy.save_utils as su
import UnitMatchPy.GUI as gui
import UnitMatchPy.assign_unique_id as aid
import UnitMatchPy.default_params as default_params
from pathlib import Path
import os


# Get default parameters, can add your own before or after!
param = default_params.get_default_param()

# snakemake.input.npy_files is a list of paths
# to get the corresponding session ids I need to parse the filenames and get the folder 2 levels up
npy_files = list(snakemake.input.npy_files)
session_ids = [Path(p).parents[2].name for p in npy_files]  # RawWaveforms -> kilosort -> session

# get unique and sorted session ids
session_ids = sorted(list(set(session_ids)))

# also get KS_dirs from the raw_waveforms paths
KS_dirs = [str(Path(p).parents[1]) for p in npy_files]  # RawWaveforms -> kilosort
KS_dirs = sorted(list(set(KS_dirs)))

param['KS_dirs'] = KS_dirs
wave_paths, unit_label_paths, channel_pos = util.paths_from_KS(KS_dirs)
param = util.get_probe_geometry(channel_pos[0], param)

# STEP 0 -- data preparation
# Read in data and select the good units and exact metadata
# if snakemake.config['unit_match']['match_good_units_only'] is False, all units will be "good" units
waveform, session_id, session_switch, within_session, good_units, param = util.load_good_waveforms(wave_paths, unit_label_paths, param, good_units_only = snakemake.config['unit_match']['match_good_units_only']) 

# param['peak_loc'] = #may need to set as a value if the peak location is NOT ~ half the spike width

# Create clus_info, contains all unit id/session related info
clus_info = {'good_units' : good_units, 'session_switch' : session_switch, 'session_id' : session_id, 
            'original_ids' : np.concatenate(good_units) }

# STEP 1
# Extract parameters from waveform
extracted_wave_properties = ov.extract_parameters(waveform, channel_pos, clus_info, param)

# STEP 2, 3, 4
# Extract metric scores
total_score, candidate_pairs, scores_to_include, predictors  = ov.extract_metric_scores(extracted_wave_properties, session_switch, within_session, param, niter  = 2)

# STEP 5
# Probability analysis
# Get prior probability of being a match
prior_match = 1 - (param['n_expected_matches'] / param['n_units']**2 ) # freedom of choose in prior prob
priors = np.array((prior_match, 1-prior_match))

# Construct distributions (kernels) for Naive Bayes Classifier
labels = candidate_pairs.astype(int)
cond = np.unique(labels)
score_vector = param['score_vector']
parameter_kernels = np.full((len(score_vector), len(scores_to_include), len(cond)), np.nan)

parameter_kernels = bf.get_parameter_kernels(scores_to_include, labels, cond, param, add_one = 1)

# Get probability of each pair of being a match
probability = bf.apply_naive_bayes(parameter_kernels, priors, predictors, param, cond)

output_prob_matrix = probability[:,1].reshape(param['n_units'],param['n_units'])

util.evaluate_output(output_prob_matrix, param, within_session, session_switch, match_threshold = 0.75)

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
gui.process_info_for_GUI(output_prob_matrix, match_threshold, scores_to_include, total_score, amplitude, spatial_decay,
                         avg_centroid, avg_waveform, avg_waveform_per_tp, wave_idx, max_site, max_site_mean, 
                         waveform, within_session, channel_pos, clus_info, param)


matches = np.argwhere(output_threshold == 1)
UIDs = aid.assign_unique_id(output_prob_matrix, param, clus_info)

# save_dir is the directory of snakemake.output.match_table
save_dir = os.path.dirname(snakemake.output.match_table)
#NOTE - change to matches to matches_curated if done manual curation with the GUI
su.save_to_output(save_dir, scores_to_include, matches # matches_curated
                  , output_prob_matrix, avg_centroid, avg_waveform, avg_waveform_per_tp, max_site,
                   total_score, output_threshold, clus_info, param, UIDs = UIDs, matches_curated = None, save_match_table = True)