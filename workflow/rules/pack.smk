import os

src_path = config['src_path']
dst_path = config['dst_path']
session  = config['sessions'][0]

# wrappers to build absolute paths to source / destination folders
abs_src = lambda file_name: os.path.join(src_path, session.split('_')[0], session, file_name)
abs_dst = lambda file_name: os.path.join(dst_path, session.split('_')[0], session, file_name)


rule pack:
    input:
        abs_src("positions.csv"),
        abs_src("events.csv"),
        abs_src("sounds.csv"),
        abs_src("%s.json" % session)
    output:
        abs_dst("meta.h5")
    script:
        "../scripts/pack.py"