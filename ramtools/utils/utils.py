import numpy as np
import os
import shutil
import tempfile
import contextlib


@contextlib.contextmanager
def temporary_cd(dir_path):
    prev_dir = os.getcwd()
    os.chdir(os.path.abspath(dir_path))
    try:
        yield
    finally:
        os.chdir(prev_dir)

def read_xvg(fname):
    data=[]
    with open(fname) as f:
        for line in f:
            # Lines with metadata or comments start with #, @
            if not line.startswith(("@","#")):
                data.append(np.array([float(s) for s in line.split()]))
    data = np.vstack(data)
    return data
