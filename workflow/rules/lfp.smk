import os
import numpy as np
import h5py
import xml.etree.ElementTree as ET


# Replace the sampling rate value for ndm lfp
rule update_lfp_rate:
    input:
        xml=n_path('{animal}', '{session}', '{session}.xml')
    output:
        xml=temp(n_path('{animal}', '{session}', '{session}.lfp.xml'))
    run:
        with open(input.xml, 'r') as f:
            filedata = f.read()

        ndm_idx = filedata.find('ndm_lfp')
        vl_idx = filedata[ndm_idx:].find('<value>')
        vr_idx = filedata[ndm_idx:].find('</value>')
        f_upd = filedata[:ndm_idx + vl_idx + 7] + str(config['lfp']['s_rate']) + filedata[ndm_idx + vr_idx:]

        # update existing and create a new temp file
        with open(input.xml, 'w') as f:
            f.write(f_upd)
        with open(output.xml, 'w') as f:
            f.write(f_upd)


# extract LFP using neurosuite
rule extract_lfp:
    input:
        xml=n_path('{animal}', '{session}', '{session}.lfp.xml'),
        dat=n_path('{animal}', '{session}', '{session}.dat')
    output:
        lfp=temp(n_path('{animal}', '{session}', '{session}.lfp'))
    params:
        session="{session}",
        animal="{animal}"
    shell:
        "cd %s; %s %s" % (
            n_path('{params.animal}', '{params.session}', ''),
            os.path.join(config['ndm_path'], "ndm_lfp"),
            '{params.session}.xml'
        )


# convert from binary to HDF5 format, shift by offset
rule lfp2hdf5:
    input:
        xml=n_path('{animal}', '{session}', '{session}.xml'),
        lfp=n_path('{animal}', '{session}', '{session}.lfp')
    output:
        lfp_h5=os.path.join(config['dst_path'], '{animal}', '{session}', 'lfp.h5')
    run:
        # determine channel numbers
        ch_num = int(ET.parse(input.xml).getroot().findall('acquisitionSystem')[0].findall('nChannels')[0].text)

        # read all LFP in one block
        block = np.fromfile(input.lfp, dtype=np.int16)

        # TODO shift by an offset!

        with h5py.File(output.lfp_h5, 'w') as f:
            lfp_ds = f.create_dataset('lfp', data=block.reshape(int(block.shape[0]/ch_num), ch_num))
            lfp_ds.attrs['headers'] = 'samples, channels'