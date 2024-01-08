# Description: This script is used to analyse the videos using the trained model.
# It is called by the snakefile and uses the config file to find the trained model.

import os
from os.path import join

os.environ["CUDA_VISIBLE_DEVICES"]=str(snakemake.config['cuda_device_to_use'])

import deeplabcut

dlc_config_path = snakemake.config['dlc_config_path'] # path of the config file for the trained model
dlc_shuffle = snakemake.config['dlc_shuffle']
dlc_trainingsetindex = snakemake.config['dlc_trainingsetindex']
 

video_path = [snakemake.input.video_path]

session = snakemake.params.session
animal = snakemake.params.animal

DLCScorer = deeplabcut.analyze_videos(dlc_config_path,video_path,shuffle=3, videotype='avi',save_as_csv=True,trainingsetindex=0,allow_growth=True, \
                                      destfolder=join(snakemake.config["dst_path"],animal,session,'dlc'))



