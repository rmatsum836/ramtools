import os
from pkg_resources import resource_filename

def get_fn(name):
    """Get the full path to one of the reference files shipped for utils.

    Parameters
    ----------
    name : str
        Name of the file to load (with respect to the reference/ folder).
    """
    fn = resource_filename('ramtools', os.path.join('tests', 'files', name))
    if not os.path.exists(fn):
        raise IOError('Sorry! {} does not exists.'.format(fn))
    return fn
