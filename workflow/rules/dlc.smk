
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

rule run_dlc_inference:
    input:
        video_path = join(config["src_path"],"{animal}","{session}","{videoname}.avi")
    output:
        DLC_h5_path = join(config["dst_path"],"{animal}","{session}","dlc","{videoname}"+DLC_scorer+".h5"),
        DLC_h5_path_meta = join(config["dst_path"],"{animal}","{session}","dlc","{videoname}"+DLC_scorer+"_meta.pickle")
    params:
        session = "{session}",
        animal = "{animal}"
    conda:
        "/mnt/nevermind.data-share/ag-grothe/AG_Pecka/envs/dlc"
    script:
        "../scripts/analyse_video.py"


ruleorder: create_labeled_processed_video > create_labeled_raw_video

rule create_labeled_raw_video:
    input:
        video_path = join(config["src_path"],"{animal}","{session}","{videoname}.avi"),
        DLC_h5_path = join(config["dst_path"],"{animal}","{session}","dlc","{videoname}"+DLC_scorer+".h5"),
        DLC_h5_path_meta = join(config["dst_path"],"{animal}","{session}","dlc","{videoname}"+DLC_scorer+"_meta.pickle")
    output:
        labeled_raw_video = join(config["dst_path"],"{animal}","{session}","{videoname}_labeled.mp4")
    params:
        session = "{session}",
        animal = "{animal}",
        videoname = "{videoname}"
    conda:
        "/mnt/nevermind.data-share/ag-grothe/AG_Pecka/envs/dlc"
    script:
        "../scripts/create_labeled_raw_video.py"

# TODO: This rule is WiP
rule create_labeled_processed_video:
    input:
        DLC_h5_path = lambda wildcards: join(config["dst_path"],"{animal}","{session}","dlc",get_videoname(wildcards,config)+DLC_scorer+".h5"),
        DLC_h5_path_meta = lambda wildcards: join(config["dst_path"],"{animal}","{session}","dlc",get_videoname(wildcards,config)+DLC_scorer+"_meta.pickle"),
        labeled_raw_video = lambda wildcards: join(config["dst_path"],"{animal}","{session}",get_videoname(wildcards,config)+"_labeled.mp4"),
        proc_video = join(config["src_path"],"{animal}","{session}","video.avi"),
    output:
        labeled_proc_video = join(config["dst_path"],"{animal}","{session}","video_labeled.mp4")
    params:
        session = "{session}",
        animal = "{animal}",
    conda:
        "/mnt/nevermind.data-share/ag-grothe/AG_Pecka/envs/dlc"
    # script:
    #     "../scripts/create_labeled_raw_video.py"
    shell:
        """
        cp {input.labeled_raw_video} {output}
        """


rule synch_dlc_100Hz:
    input:
        DLC_h5_path = lambda wildcards: join(config["dst_path"],"{animal}","{session}","dlc",get_videoname(wildcards,config)+DLC_scorer+".h5"),
        events = join(config["src_path"],"{animal}","{session}","events.csv"),
        session_cfg = join(config["src_path"],"{animal}","{session}","{session}.json"),
        timestamps_raw = lambda wildcards: join(config["src_path"],"{animal}","{session}",get_videoname(wildcards,config)+".csv"),
    output:
        DLC_h5_path_100Hz = join(config["dst_path"],"{animal}","{session}","DLC_100Hz.h5"),
    conda:
        "/mnt/nevermind.data-share/ag-grothe/AG_Pecka/envs/dlc"
    script:
        "../scripts/synch_DLC_100Hz.py"

rule create_video_timestamps:
    input:
        video_path = lambda wildcards: join(config["src_path"],"{animal}","{session}","video.avi"),
    output:
        timestamps = join(config["src_path"],"{animal}","{session}","video.csv"),
    run:
        import cv2
        import numpy as np

        video = cv2.VideoCapture(input.video_path)

        # Get the frames per second
        fps = video.get(cv2.CAP_PROP_FPS)

        # Get the total number of frames
        frame_count = int(video.get(cv2.CAP_PROP_FRAME_COUNT))

        video.release()

        # Create a np array with the timestamps
        timestamps = np.linspace(0, frame_count/fps, frame_count)

        # Save the timestamps to a csv file
        np.savetxt(output.timestamps, timestamps, delimiter="\n")
