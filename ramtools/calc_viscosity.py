import numpy as np
import mdtraj as md
import scipy
import unyt as u
import os
from ramtools.utils.utils import read_xvg, temporary_cd


def get_energies(energy_file, volume):
    if not os.path.exists('energy.xvg'):
        cmd = f'echo {volume} | gmx energy -f {energy_file} -vis'
        os.system(cmd)

def calc_green_kubo(trj_file, top_file, energy_file, temp=300):
    """ Calculate Green-Kubo Viscosity
    """

    temp = temp * u.K

    trj = md.load(trj_file, top=top_file)
    volume = np.mean(trj.unitcell_volumes) * u.nm**3

    get_energies(energy_file, volume)
    data = read_xvg('energy.xvg')[:50000]
    print("Data correctly loaded")
    time = (data[:,0] * u.ps).to(u.s)
    xy = (data[:,2] * u.bar).to(u.pascal) 
    xz = (data[:,3] * u.bar).to(u.pascal)
    yz = (data[:,7] * u.bar).to(u.pascal)
   
    pressures = list()
    for i in range(len(xy)):
        xy_mult = xy[i] * xy[0]
        xz_mult = xz[i] * xz[0]
        yz_mult = yz[i] * yz[0]
        pressures.append(np.mean([xy_mult, xz_mult, yz_mult])*u.Pa**2)

    #xy = [np.mean(xy[i] * xy[0]) for i in range(len(xy))]
    joules = (u.kg * u.m**2) / u.s**2

    kb = 1.38e-23 * joules

    volume = volume.to(u.m**3)
    coefficient = volume / (kb * temp)

    #pascal = u.kg / (u.m * u.s**2)
    integral = [np.trapz(pressures[:i], time[:i]) for i in range(len(pressures))]
    #integral = [pressures[i]*time[i] for i in range(len(pressures))]
    #integral *= pascal**2 * u.s
    return time, [val * coefficient  for val in integral]

def calc_multiple_green_kubo(trj_file, top_file, energy_file, n_runs, temp=300):
    """ Calculate Green-Kubo Viscosity From Multiple Trajectories
    """
    integral_list = list()
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
        integral_list.append(integral)


    # Take average of running integrals
    average = np.mean(integral_list, axis=0)  
    return time, [val * coefficient  for val in average]

def calc_multiple_einstein(trj_file, top_file, energy_file, n_runs):
    """ Calculate Einstein Viscosity From Multiple Trajectories
    """
    rho_list = list()
    for i in range(1, n_runs):
        dirname = 'run_{}'.format(i)
        if not os.path.isdir(dirname):
            raise OSError("Directory doesn't exist")
        
        with temporary_cd(dirname):
            trj = md.load(trj_file, top=top_file)
            volume = np.mean(trj.unitcell_volumes)
            get_energies(energy_file, volume)
            data = read_xvg('evisco.xvg')
            print("Data correctly loaded")
            time = data[:,0]
            rhos = np.mean([data[:,1], data[:,2], data[:,3], data[:,4]], axis=0)

        rho_list.append(rhos)


    # Take average of running integrals
    average = np.mean(rho_list, axis=0)  
    return time, average
