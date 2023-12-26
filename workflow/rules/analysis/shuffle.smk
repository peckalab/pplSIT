

# main function to create all inputs and wildcards
def get_psth_micro_inputs(w):
    units_file = os.path.join(config['dst_path'], w.animal, w.session, 'units.h5')
    with h5py.File(units_file, 'r') as f:
        unit_names = [name for name in f]  
    return expand(os.path.join(config['dst_path'], w.animal, w.session, 'analysis', '{unit_name}_psth_shuffle_micro_silence.npy'), unit_name=unit_names)


rule shuffle_micro:
    input:
        meta=os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        units=os.path.join(config['dst_path'], '{animal}', '{session}', 'units.h5')
    output:
        dst_unit_file=temp(os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', '{unit_name}_psth_shuffle_micro_silence.npy'))
    script:
        "../../scripts/analysis/shuffle.py"


rule shuffle_micro_pack:
    input:
        get_psth_micro_inputs
    output:
        dst_file=os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'psth_shuffle_micro_silence.h5')
    run:
        for input_file in input:
            unit_name = os.path.basename(input_file).split('_')[0]
            with h5py.File(output.dst_file, 'a') as f:
                if unit_name in f:
                    del f[unit_name]
                f.create_dataset(unit_name, data=np.load(input_file))