import numpy as np
import mdtraj as md
import scipy
import os
from ramtools.utils.utils import read_xvg

def get_energies(energy_file, volume):
    if not os.path.exists('energy.xvg'):
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
    time = data[:,0]
    xy = data[:,2]
    xz = data[:,3]
    yz = data[:,7]
    #xy *= 100000
    #xy *= 100000
    #xy *= 100000
   
    pressures = list()
    for i in range(len(xy)):
        xy_mult = xy[i] * xy[0]
        xz_mult = xz[i] * xz[0]
        yz_mult = yz[i] * yz[0]
        pressures.append(np.mean([xy_mult, xz_mult, yz_mult]))

    #xy = [np.mean(xy[i] * xy[0]) for i in range(len(xy))]

    kb = 1.38e-23

    volume *= 1e-27
    coefficient = volume / (kb * temp)


    #integral = [np.trapz(xy[:i], time[:i]) for i in range(len(xy))]
    integral = [np.trapz(pressures[:i], time[:i]) for i in range(len(pressures))]
  
    return time, [val * coefficient  for val in integral]
