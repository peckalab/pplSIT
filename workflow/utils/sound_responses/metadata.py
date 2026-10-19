"""Stable method and reference identifiers shared by both stages."""
METHOD = 'amaro_07062021_bayesian_blocks'
METHOD_VERSION = 1
CACHE_SCHEMA_VERSION = 1
CACHE_ARTIFACT = 'sound_responses_bayesian_blocks_segmentations'

REFERENCE = {
    'file': 'code_07062021/neurofunctions/general.py',
    'functions': ['determine_latency_bayesianblocks', 'determine_offset'],
    'sha256': 'c96d5c851fed1723d950eab61f86b93d0e9d8974f078f03a6e9a220b4df02735',
    'backend_reference': 'skhep.modeling.bayesian_blocks (scikit-hep 0.3.1 comparison)',
}
