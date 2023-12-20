import os
import json


rule extract_aeps:
    input:
        manual=os.path.join(config['src_path'], '{animal}', '{session}', 'manual.json'),
        lfp=os.path.join(config['dst_path'], '{animal}', '{session}', 'lfp.h5'),
        meta=os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
    output:
        aeps=os.path.join(config['dst_path'], '{animal}', '{session}', 'AEPs.h5')
    script:
        "../scripts/aeps.py"
        