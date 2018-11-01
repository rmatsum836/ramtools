import numpy as np
import os
import sys
import mdtraj as md
from mtools.gromacs.gromacs import make_comtrj
import matplotlib as mpl
import matplotlib.pyplot as plt

def calc_charge_density(gro_file, trj_file, top_file, bin_width, area,
        dim, box_range, maxs, mins):

    trj = md.load(trj_file, top=gro_file)
    com_trj = make_comtrj(trj)
    resnames = np.unique([x.name for x in
               com_trj.topology.residues])
    bin_list = list()
    q_list = list()
