print("DLC_h5_path: ", snakemake.input.DLC_h5_path)
print("moseq_csv_path: ", snakemake.output.moseq_csv_path)

import os, glob
from os.path import join

os.environ["CUDA_VISIBLE_DEVICES"]=str(snakemake.config['cuda_device_to_use'])

os.environ['TF_FORCE_GPU_ALLOW_GROWTH'] = 'true'
os.environ['XLA_PYTHON_CLIENT_MEM_FRACTION'] = '0.75'

import keypoint_moseq as kpms

moseq_model_path = snakemake.config['moseq_model_path']

project_dir = os.path.dirname(moseq_model_path)
config = lambda: kpms.load_config(project_dir)

name = os.path.basename(moseq_model_path)

# # load the most recent model checkpoint and pca object
checkpoint = kpms.load_checkpoint(project_dir, name, path = join(moseq_model_path,'checkpoint.p'))
pca = kpms.load_pca(project_dir)

# # load new data (e.g. from deeplabcut)
new_data = snakemake.input.DLC_h5_path # can be a file, a directory, or a list of files
dlc_h5_name = os.path.basename(new_data)
dlc_h5_name_without_extension = os.path.splitext(dlc_h5_name)[0]
coordinates, confidences, bodyparts = kpms.load_keypoints(new_data, 'deeplabcut')


results = kpms.apply_model(
    coordinates=coordinates, 
    confidences=confidences, 
    project_dir=project_dir, 
    name=name, 
    pca=kpms.load_pca(project_dir),
    params=checkpoint['params'],
    hypparams=checkpoint['hypparams'],
    **config())

# optionally rerun `save_results_as_csv` to export the new results
kpms.save_results_as_csv(project_dir=project_dir,name=name,save_dir=join(snakemake.config["dst_path"],snakemake.params.animal,snakemake.params.session,'moseq'))

# change name of the file to moseq.csv
os.rename(join(snakemake.config["dst_path"],snakemake.params.animal,snakemake.params.session,'moseq',dlc_h5_name_without_extension+'.csv'),\
          snakemake.output.moseq_csv_path)

# Get a list of all csv files in the directory
csv_files = glob.glob(join(snakemake.config["dst_path"], snakemake.params.animal, \
                           snakemake.params.session, 'moseq', '*.csv'))

# Remove all csv files except the one we want to keep
for file in csv_files:
    if file != snakemake.output.moseq_csv_path:
        os.remove(file)
