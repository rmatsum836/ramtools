import numpy as np
import mdtraj as md
import scipy
import os
from utils.utils import read_xvg
from scipy.integrate import simps

def get_energies(energy_file, volume):
    cmd = f'echo {volume} | gmx energy -f {energy_file} -vis'
    os.system(cmd)

def calc_green_kubo(trj_file, top_file, energy_file, temp=300):
    """ Calculate Green-Kubo Viscosity
    """

    trj = md.load(trj_file, top=top_file)
    volume = np.mean(trj.unitcell_volumes)

    get_energies(energy_file, volume)
    data = read_xvg('energy.xvg')
    print("Data correctly loaded")
    xy = data[:,2]
    xz = data[:,3]
    yz = data[:,7]
   
    #pressures = np.empty(shape=(len(xy), 3))
    #for i in pressures:
    #    pressures[i] = np.array(xy[i], xz[i], yz[i])
    #xy_mean np.mean([xy[0]*i for i in xy])

    kb = 1.38e-23

    volume *= 1e-27
    pressures *= 100000
    coefficient = volume / (kb * T)

    integral = simps(xy[t]*xy[0], t)

    return coefficient * integral
