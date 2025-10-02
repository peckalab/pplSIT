import hashlib


def pval2text(p_val):
    if p_val > 0.05:
        return 'n.s.'
    elif p_val > 0.01:
        return '*'
    elif p_val > 0.001:
        return '**'
    elif p_val > 0.0001:
        return '***'
    elif p_val > 0.00001:
        return '****'
    else:
        return '*****'
    

def hex_hash(text, length=10):
    return hashlib.sha1(text.encode()).hexdigest()[:length]