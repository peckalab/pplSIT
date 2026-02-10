from copy import deepcopy

def recursive_merge(base, override):
    """Recursively merge override dict into base dict"""
    result = deepcopy(base)
    for key, value in override.items():
        if key in result and isinstance(result[key], dict) and isinstance(value, dict):
            result[key] = recursive_merge(result[key], value)
        else:
            result[key] = value
    return result