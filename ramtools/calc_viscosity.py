import numpy as np
import mdtraj as md
import scipy
import os
from ramtools.utils.utils import read_xvg, temporary_cd

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

def calc_multiple_green_kubo(trj_file, top_file, energy_file, n_runs, temp=300):
    """ Calculate Green-Kubo Viscosity From Multiple Trajectories
    """
    xy_list = list()
    xz_list = list()
    yz_list = list()
    for i in range(n_runs):
        dirname = 'run_{}'.format(i)
        if not os.path.isdir(dirname):
            raise OSError("Directory doesn't exist")
        
        with temporary_cd(dirname):
            trj = md.load(trj_file, top=top_file)
            volume = np.mean(trj.unitcell_volumes)
            get_energies(energy_file, volume)
            data = read_xvg('energy.xvg')
            print("Data correctly loaded")
            time = data[:,0]
            xy = data[:,2]
            xz = data[:,3]
            yz = data[:,7]
            xy_list.append(xy)
            xz_list.append(xz)
            yz_list.append(yz)
        
    xy_mean = np.mean(xy_list, axis=0)
    xz_mean = np.mean(xz_list, axis=0)
    yz_mean = np.mean(yz_list, axis=0)
     
    pressures = list()
    for i in range(len(xy)):
        xy_mult = xy_mean[i] * xy_mean[0]
        xz_mult = xz_mean[i] * xz_mean[0]
        yz_mult = yz_mean[i] * yz_mean[0]
        pressures.append(np.mean([xy_mult, xz_mult, yz_mult]))

    #xy = [np.mean(xy[i] * xy[0]) for i in range(len(xy))]

    kb = 1.38e-23

    volume *= 1e-27
    coefficient = volume / (kb * temp)


    #integral = [np.trapz(xy[:i], time[:i]) for i in range(len(xy))]
    integral = [np.trapz(pressures[:i], time[:i]) for i in range(len(pressures))]
  
    return time, [val * coefficient  for val in integral]
