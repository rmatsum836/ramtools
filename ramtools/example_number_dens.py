import numpy as np
import os
import sys
import mdtraj as md
from mtools.gromacs.gromacs import make_comtrj
import matplotlib as mpl
import matplotlib.pyplot as plt

gro_file = 'sample.gro'
trj_file = 'sample.trr'
top_file = 'init.top'
box_range = [.92, 2.119] #2.119
mins = [12.75,0,0]
maxs = [15.583,0,0]
dim = 1
bin_width = 0.01
area = 2.973 * 2.702
rho, bins = calc_number_density(gro_file, trj_file, top_file,
        bin_width, area, dim, box_range, maxs, mins)
plt.figure()
plt.xlim(1,1.8)
plt.plot(bins[0],rho[0])
plt.plot(bins[1], rho[1])
#plt.plot(bins[2], rho[2])
plt.savefig('number-density.pdf')
