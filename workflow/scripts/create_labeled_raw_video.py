import os

from os.path import join

import shutil

os.environ["CUDA_VISIBLE_DEVICES"]=str(snakemake.config['cuda_device_to_use'])

import deeplabcut

def get_DLCScorer(dlc_h5_file_path):
        """
        Retrieves the DLC scorer from the specified DLC H5 file.

        Parameters:
        dlc_h5_file_path (str): The file path to the DLC H5 file.

        Returns:
        str: The name of the DLC scorer.

        """
        import pandas as pd
        df_dlc = pd.read_hdf(dlc_h5_file_path)
        DLCScorer = df_dlc.columns[0][0]
        return DLCScorer


dlc_config_path = snakemake.config['dlc_config_path'] # path of the config file for the trained model
dlc_shuffle = snakemake.config['dlc_shuffle']
dlc_trainingsetindex = snakemake.config['dlc_trainingsetindex']


video_path = [snakemake.input.video_path]

session = snakemake.params.session
animal = snakemake.params.animal
destfolder = join(snakemake.config["dst_path"],snakemake.params.animal,snakemake.params.session,'dlc')


# # TODO: this is a very ugly workaround to format the h5 filename with the DLC scorer name like 
# # create_labeled_video expects it
DLCScorer = get_DLCScorer(snakemake.input.DLC_h5_path)
reformatted_h5_filepath = join(destfolder,"raw"+DLCScorer+".h5")
# shutil.copyfile(snakemake.input.DLC_h5_path, reformatted_h5_filepath)

deeplabcut.create_labeled_video(dlc_config_path,video_path,destfolder=destfolder, \
                                shuffle=dlc_shuffle, videotype='avi', filtered=False,\
                                trainingsetindex=dlc_trainingsetindex,save_frames=False,keypoints_only=False)

# os.remove(reformatted_h5_filepath)

os.rename(join(snakemake.config["dst_path"],animal,session,'dlc',snakemake.params.videoname+DLCScorer+"_labeled.mp4"), \
          snakemake.output.labeled_raw_video) # rename and move the output file to the desired location