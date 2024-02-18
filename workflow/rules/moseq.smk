import os
from os.path import join

# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(parent_dir)
sys.path.append(os.path.join(parent_dir, 'utils'))
from utils.dlc import get_scorer_name, read_config

DLC_scorer,_ = get_scorer_name(read_config(config["dlc_config_path"]),config["dlc_shuffle"],config["dlc_trainingsetindex"])


def get_videoname(wildcards, config):
    """
    Returns the video name based on the existence of a raw.avi file.

    Args:
        wildcards (object): The wildcards object containing the animal and session information.
        config (dict): The configuration dictionary containing the source path.

    Returns:
        str: The video name ("raw" or "video").
    """

    raw_path = join(config["src_path"], wildcards.animal, wildcards.session, "raw.avi")
    return "raw" if os.path.exists(raw_path) else "video"


rule apply_moseq:
    input:
        DLC_h5_path = lambda wildcards: join(config["dst_path"],"{animal}","{session}","dlc",get_videoname(wildcards,config)+DLC_scorer+".h5"),
    output:
        moseq_csv_path = join(config["dst_path"],"{animal}","{session}","moseq","moseq.csv"),
    params:
        session = "{session}",
        animal = "{animal}"
    conda:
        "/mnt/nevermind.data-share/ag-grothe/AG_Pecka/envs/kpms_v0_1_5"
    script:
        "../scripts/apply_moseq_new_data.py"


# TODO: add rule to synch moseq data to 100Hz
# rule synch_moseq_100Hz:
#     input:
#         moseq_csv = join(config["dst_path"],"{animal}","{session}","moseq","moseq.csv"),
#         events = join(config["src_path"],"{animal}","{session}","events.csv"),
#         session_cfg = join(config["src_path"],"{animal}","{session}","{session}.json"),
#         timestamps_raw = lambda wildcards: join(config["src_path"],"{animal}","{session}",get_videoname(wildcards,config)+".csv"),
#     output:
#         DLC_h5_path_100Hz = join(config["dst_path"],"{animal}","{session}","DLC_100Hz.h5"),
#     conda:
#         "/mnt/nevermind.data-share/ag-grothe/AG_Pecka/envs/kpms_v0_1_5"
#     script:
#         "../scripts/synch_moseq_100Hz.py"