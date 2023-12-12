import os

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