import xml.etree.ElementTree as ET
import os, sys

# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(parent_dir)
sys.path.append(os.path.join(parent_dir, 'utils'))
from utils.neurosuite import XMLHero


#def build_ndm_targets(f_type, animal, session):
#    return [n_path(animal, session, '%s.%s.%s' % (session, f_type, el)) for el in [1, 2, 3, 4, 5, 6, 7, 8]]

# main function to create all inputs and wildcards
def get_clu_inputs(w):
    xml_path = os.path.join(config['src_path'], w.animal, w.session, w.session + '.xml')
    electrodes = XMLHero(xml_path).get_electrodes()
    return expand(n_path(w.animal, w.session, w.session + '.clu_copy.' + '{electrode}'), electrode=electrodes)

def get_fet_param(w):
    xml_path = os.path.join(config['src_path'], w.animal, w.session, w.session + '.xml')
    return XMLHero(xml_path).get_fet_string(w.electrode)

# -------- rules ----------------------


# high-pass filtering
rule hipass:
    input:
        xml=n_path('{animal}', '{session}', '{session}.xml'),
        dat=n_path('{animal}', '{session}', '{session}.dat')
    output:
        temp(n_path('{animal}', '{session}', '{session}.fil'))
    params:
        session="{session}",
        animal="{animal}"
    shell:
        "cd %s; %s %s" % (
            n_path('{params.animal}', '{params.session}', ''),
            os.path.join(config['ndm_path'], "ndm_hipass"),
            '{params.session}.xml'
        )

# extract spikes
# to run the rule only once, we require only one dummy 
# output file 'spk.ready' as output instead of all .spk / .res
# files. No better way to run it with snakemake found
rule extractspikes:
    input:
        xml=n_path('{animal}', '{session}', '{session}.xml'),
        fil=n_path('{animal}', '{session}', '{session}.fil')
    output:
        temp(n_path('{animal}', '{session}', 'spk.ready'))
        #expand(n_path('{animal}', '{session}', '{session}.spk.{el}'), el=ACTUALS['electrodes']),
        #expand(n_path('{animal}', '{session}', '{session}.res.{el}'), el=ACTUALS['electrodes'])
        #build_ndm_targets('spk', '{animal}', '{session}'),
        #build_ndm_targets('res', '{animal}', '{session}')
        #n_path('{animal}', '{session}', '{session}.spk.{electrode}'),  # this works, but makes the 
        #n_path('{animal}', '{session}', '{session}.res.{electrode}')   # rule be executed multiple times
    params:
        session="{session}",
        animal="{animal}"
    shell:
        "cd %s; %s %s; touch %s" % (
            n_path('{params.animal}', '{params.session}', ''),
            os.path.join(config['ndm_path'], "ndm_extractspikes"),
            '{params.session}.xml',
            n_path('{params.animal}', '{params.session}', 'spk.ready')
        )

# compute PCA
# to run the rule only once, we require only one dummy 
# output file 'fet.ready' as output instead of all .fet
# files. No better way to run it with snakemake found
rule PCA:
    input:
        #n_path('{animal}', '{session}', '{session}.spk.{electrode}')
        n_path('{animal}', '{session}', 'spk.ready')
    output:
        #n_path('{animal}', '{session}', '{session}.fet.{electrode}')
        temp(n_path('{animal}', '{session}', 'fet.ready'))
    params:
        session="{session}",
        animal="{animal}"
    shell:
        "cd %s; %s %s; touch %s" % (
            n_path('{params.animal}', '{params.session}', ''),
            os.path.join(config['ndm_path'], "ndm_pca"),
            '{params.session}.xml',
            n_path('{params.animal}', '{params.session}', 'fet.ready')
        )


# sort clusters with KlustaKwik
rule kkwik:
    input:
        #n_path('{animal}', '{session}', '{session}.fet.{electrode}'),
        #n_path('{animal}', '{session}', '{session}.res.{electrode}')
        n_path('{animal}', '{session}', 'fet.ready')
    output:
        n_path('{animal}', '{session}', '{session}.clu.{electrode}')
    params:
        session="{session}",
        animal="{animal}",
        electrode="{electrode}",
        fet_string=lambda w, output: get_fet_param(w)
    shell:
        "cd %s; %s %s %s %s %s" % (
            n_path('{params.animal}', '{params.session}', ''),
            os.path.join(config['kkwik_path'], "KlustaKwik"),
            '{params.session}',
            '{params.electrode}',
            config['kwik_args'],
            '{params.fet_string}'
        )


# backup .clu files
rule clu_backup:
    input:
        n_path('{animal}', '{session}', '{session}.clu.{electrode}')
    output:
        n_path('{animal}', '{session}', '{session}.clu_copy.{electrode}')
    shell:
        "cp {input} {output}"


 # dump units to HDF5
rule dump_units:
    input:
        xml=n_path('{animal}', '{session}', '{session}.xml'),
        meta=os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        clu=get_clu_inputs
    output:
        os.path.join(config['dst_path'], '{animal}', '{session}', 'units.h5')
    script:
        "../scripts/units.py"
