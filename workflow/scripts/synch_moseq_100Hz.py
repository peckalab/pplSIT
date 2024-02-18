from os.path import join
import pandas as pd
import numpy as np
import h5py
import json

def px_to_meters(cfg, x, y):  # convert pixels to meters
    cfg_pos = cfg['position']
    pixel_size = cfg_pos['floor_r_in_meters'] / float(cfg_pos['floor_radius'])
    x_m = (cfg_pos['arena_x'] - x).astype(float) * pixel_size * (-1 if cfg_pos['flip_x'] else 1)
    y_m = (cfg_pos['arena_y'] - y).astype(float) * pixel_size * (-1 if cfg_pos['flip_y'] else 1)
    return x_m, y_m

with open(snakemake.input.session_cfg, 'r') as f:
    session_cfg = json.load(f)

df_moseq_30Hz = pd.read_csv(snakemake.input.moseq_csv_path)

df_moseq_100Hz = pd.DataFrame(columns=['time']+[column for column in df_moseq_30Hz.columns])

events = np.loadtxt(snakemake.input.events, delimiter=',', skiprows=1)

time_freq = 100  # at 100Hz

s_start, s_end = events[:, 0][0], events[:, 0][-1]
events[:, 0] = events[:, 0] - s_start
times = np.arange(0, events[:, 0][-1], (1./ time_freq))
timestamps_raw = np.loadtxt(snakemake.input.timestamps_raw)
timestamps_raw = timestamps_raw - s_start

# Get fps of video
with open(snakemake.input.session_cfg) as json_file:
    parameters = json.load(json_file)
fps_video = parameters['video']['fps']

# Find the indices of the timestamps that are closest to the times
indices = np.searchsorted(timestamps_raw, times)

df_moseq_100Hz['time'] = times
for column in df_moseq_30Hz.columns:
    if "syllables" not in column:
        df_moseq_100Hz[column] = np.interp(times,timestamps_raw,df_moseq_30Hz[column])
    else:
        # clip indices to be within the range of the original data (between 0 and len(df_moseq_30Hz[column]) - 1)
        indices = np.clip(indices, 0, len(df_moseq_30Hz[column]) - 1)
        df_moseq_100Hz[column] = df_moseq_30Hz[column].iloc[indices].values


array_moseq_100Hz = df_moseq_100Hz.to_numpy()

for col_idx, column in enumerate(df_moseq_100Hz.columns):
    if ' x' in column:
        array_moseq_100Hz[:,col_idx],array_moseq_100Hz[:,col_idx+1] = px_to_meters(session_cfg,array_moseq_100Hz[:,col_idx],array_moseq_100Hz[:,col_idx+1])

with h5py.File(snakemake.output.moseq_h5_path_100Hz, 'w') as f:
    df_moseq_dataset = f.create_dataset('moseq',data=array_moseq_100Hz)
    df_moseq_dataset.attrs['headers'] = np.array(df_moseq_100Hz.columns.values)