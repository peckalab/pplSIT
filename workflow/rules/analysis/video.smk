import os
import h5py


def get_video_inputs(w):
    meta_file = os.path.join(config['dst_path'], w.animal, w.session, 'meta.h5')
    with h5py.File(meta_file, 'r') as f:
        tl = np.array(f['processed']['timeline'])

    step = config['video']['duration'] - config['video']['overlap']
    t_starts = [step * i for i in range(round(tl[:, 0][-1] / step))]
    t_ends   = [t + config['video']['duration'] for t in t_starts[:-1]]
    t_ends.append(int(tl[:, 0][-1]))

    time_ranges = ["%s_%s" % (t1, t2) for t1, t2 in zip(t_starts, t_ends)]
    return expand(os.path.join(config['dst_path'], w.animal, w.session, 'video', '{time_range}' + '.mp4'), time_range=time_ranges)


rule crop_original:
    input:
        f=os.path.join(config['src_path'], '{animal}', '{session}', 'video.avi'),
    output:
        #f=temp(os.path.join(config['dst_path'], '{animal}', '{session}', 'video', 'video_cropped.avi'))
        f=os.path.join(config['dst_path'], '{animal}', '{session}', 'video', 'video_cropped.avi')
    shell:
        'ffmpeg -i {input.f} -vf "crop=%s:%s:0:0" -c:v libx264 -crf 17 -c:a copy {output.f}' % (
            config['video']['crop_width'], config['video']['crop_height']
        )
        

rule scale_cropped:
    input:
        f=os.path.join(config['dst_path'], '{animal}', '{session}', 'video', 'video_cropped.avi')
    output:
        f=os.path.join(config['dst_path'], '{animal}', '{session}', 'video', 'video_scaled.avi')
    shell:
        'ffmpeg -i {input.f} -vf scale=-2:%s -c:v libx264 -crf 17 {output.f}' % config['video']['scale_height']


rule compile_video:
    input:
        meta=os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        units=os.path.join(config['dst_path'], '{animal}', '{session}', 'units.h5'),
        labels=os.path.join(config['dst_path'], '{animal}', '{session}', 'video', 'sorted_unit_labels.json'),
        video=os.path.join(config['dst_path'], '{animal}', '{session}', 'video', 'video_scaled.avi')
    output:
        result=os.path.join(config['dst_path'], '{animal}', '{session}', 'video', '{time_range}.mp4')
    script:
        "../../scripts/analysis/video.py"

 
rule build_all_videos:
    input:
        meta=os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        inputs=get_video_inputs
    output:
        ready=os.path.join(config['dst_path'], '{animal}', '{session}', 'video', 'videos.ready')
    run:
        # dummy step to request all videos
        with open(output.ready, 'w') as f:
            f.write('')