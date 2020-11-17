import numpy as np

def read_xvg(fname):
    """Read an xvg file

    Parameters
    ----------
    fname : str
        filename of xvg file

    Returns
    -------
    data : np.array
        NumPy array of data
    """
    data=[]
    with open(fname) as f:
        for line in f:
            # Lines with metadata or comments start with #, @
            if not line.startswith(("@","#")):
                data.append(np.array([float(s) for s in line.split()]))
    data = np.vstack(data)
    return data
