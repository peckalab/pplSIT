import os
import json

# TODO: move this to pack - make first 'meta_nosync.h5'
# and then just 'meta.h5'?

# that's the rule to synchronize events to actual ephys
rule calibrate_event_sync:
    input:
        manual=os.path.join(config['src_path'], '{animal}', '{session}', 'manual.json'),
        meta=os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5')
    output:
        ready=temp(os.path.join(config['dst_path'], '{animal}', '{session}', 'manual.ready'))
    run:
        # original events
        with h5py.File(input.meta, 'r') as f:
            sound_events = np.array(f['processed']['sound_events'])

        # get manually detected offset
        with open(input.manual) as json_file:
            offset = json.load(json_file)['ephys']['offset']

        # update events and overwrite in meta file
        sound_events[:, 0] = sound_events[:, 0] + offset/1000.
        with h5py.File(input.meta, 'a') as f:
            del f['processed']['sound_events']
            f['processed'].create_dataset('sound_events', data=sound_events)

        # create dummy file to complete
        with open(output.ready, 'w') as f:
            f.write('foo')


rule extract_aeps:
    input:
        manual=os.path.join(config['src_path'], '{animal}', '{session}', 'manual.json'),
        lfp=os.path.join(config['dst_path'], '{animal}', '{session}', 'lfp.h5'),
        meta=os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        ready=os.path.join(config['dst_path'], '{animal}', '{session}', 'manual.ready')
    output:
        aeps=os.path.join(config['dst_path'], '{animal}', '{session}', 'AEPs.h5')
    script:
        "../scripts/aeps.py"
        