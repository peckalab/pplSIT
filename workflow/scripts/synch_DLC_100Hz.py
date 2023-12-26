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

df_DLC_30Hz = pd.read_hdf(snakemake.input.DLC_h5_path)

df_DLC_100Hz = pd.DataFrame(columns=['time']+[column for column in df_DLC_30Hz.columns])

events = np.loadtxt(snakemake.input.events, delimiter=',', skiprows=1)

time_freq = 100  # at 100Hz

events[:, 0] = events[:, 0] - events[:, 0][0]
s_start, s_end = events[:, 0][0], events[:, 0][-1]
times = np.arange(0, s_end, (1./ time_freq))
timestamps_raw = np.loadtxt(snakemake.input.timestamps_raw)

df_DLC_100Hz['time'] = times
for column in df_DLC_30Hz.columns:
    df_DLC_100Hz[column] = np.interp(times,timestamps_raw,df_DLC_30Hz[column])


df_DLC_100Hz.columns = ['time']+[column[1]+'_'+column[2] for column in df_DLC_100Hz.columns if column != 'time']

array_DLC_100Hz = df_DLC_100Hz.to_numpy()

for col_idx, column in enumerate(df_DLC_100Hz.columns):
    if '_x' in column:
        array_DLC_100Hz[:,col_idx],array_DLC_100Hz[:,col_idx+1] = px_to_meters(session_cfg,array_DLC_100Hz[:,col_idx],array_DLC_100Hz[:,col_idx+1])

with h5py.File(snakemake.output.DLC_h5_path_100Hz, 'w') as f:
    df_dlc_dataset = f.create_dataset('df_DLC',data=array_DLC_100Hz)
    df_dlc_dataset.attrs['headers'] = np.array(df_DLC_100Hz.columns.values)