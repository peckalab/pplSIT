
rule plot_spiking_maps:
    input:
        meta=os.path.join(config['dst_path'], '{animal}', '{session}', 'meta.h5'),
        units=os.path.join(config['dst_path'], '{animal}', '{session}', 'units.h5')
    output:
        os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'spiking_maps.png')
    script:
        "../../scripts/analysis/spiking_maps.py"


rule plot_firing_rate_maps:
    input:
        units=os.path.join(config['dst_path'], '{animal}', '{session}', 'units.h5')
    output:
        os.path.join(config['dst_path'], '{animal}', '{session}', 'analysis', 'fr_maps.pdf')
    script:
        "../../scripts/analysis/fr_maps.py"

#rule plot_AEP_component_maps:
