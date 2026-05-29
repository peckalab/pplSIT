import os


rule strf_passive:
    input:
        meta=os.path.join(config["dst_path"], "{animal}", "{session}", "meta.h5"),
        units=os.path.join(config["dst_path"], "{animal}", "{session}", "units.h5")
    output:
        h5=os.path.join(config["dst_path"], "{animal}", "{session}", "analysis", "strf_passive.h5"),
        pdf=os.path.join(config["dst_path"], "{animal}", "{session}", "analysis", "strf_passive.pdf")
    script:
        "../../scripts/analysis/strf_passive.py"
