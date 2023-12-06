import os

src_path = config['src_path']
dst_path = config['dst_path']
sessions = config['sessions']

def build_file_path(where, session, file_name):
    return os.path.join(where, session.split('_')[0], session, file_name)

# raw session files
raw_files = [
    build_file_path(src_path, s, '%s.csv' % ds_name)\
        for ds_name in ['positions', 'events', 'sounds', 'islands']\
        for s in sessions
]

# raw config files
cfg_files = [build_file_path(src_path, s, '%s.json' % s) for s in sessions]


# output files
out_files = [build_file_path(dst_path, s, 'meta.h5') for s in sessions]


rule pack:
    input:
        raw_files + cfg_files
    output:
        out_files
    script:
        "../scripts/pack.py"