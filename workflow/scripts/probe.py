import os, sys
import json
import shutil

# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.probe import oe2kilosort


# try to read probe settings from Open Ephys XML file
kilosort_probe = None
session_path  = os.path.dirname(snakemake.output[0])
for dirpath, dirnames, filenames in os.walk(session_path):
    for filename in [f for f in filenames if f == 'settings.xml']:
        oe_xml_path = os.path.join(dirpath, filename)

        try:
            kilosort_probe = oe2kilosort(oe_xml_path)  # this is a dict in kilosort probe format
        except ValueError:
            pass
        break

if kilosort_probe is None:  # probe description not found, copy template
    print('USING TEMPLATE FOR KILOSORT PROBE DESCRIPTION')
    shutil.copy(snakemake.input[0], snakemake.output[0])
else:
    with open(snakemake.output[0], 'w') as f:
        f.write(json.dumps(kilosort_probe, indent=2))