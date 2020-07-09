import numpy as np
import mdtraj as md
import unyt as u

from ramtools.utils.read_files import read_xvg


def calc_conductivity(N, V, D_cat, D_an, q=1, T=300):
    """ Calculate Nernst-Einstein Conductivity

    Parameters
    ----------
    N : int
        Number of ions
    V : float
        Volume of simulation box
    D_cat : float
        Diffusivity of cation in m^2/s
    D_an : float
        Diffusivity of anion in m^2/s
    q : float, default=1
        Charge of ions in element charge units
    T : float, default=300
        Temperature of system in Kelvin

    Returns
    -------
    cond : unyt.array
        Nernst-Einstein conductivity

    """
    D_cat *= u.m**2 / u.s
    D_an *= u.m**2 / u.s
    kT = T * 1.3806488e-23 * u.joule
    q *= u.elementary_charge
    q = q.to('Coulomb')

    cond = N / (V*kT) * q ** 2 * (D_cat + D_an)

    return cond


def calc_hfshear(energy_file, trj, temperature):
    """ Calculate High-Frequency shear modulus of an MDTraj trajectory

    Parameters
    ----------
    energy_file : str
        GROMACS .edr file
    trj : str
        MDTraj trajectory
    temperatrue : flt
        Temperature in Kelvin

    Returns
    -------
    shear_bar : unyt.array
        average shear modulus
    shear_std : unyt_array
        stadndard deviation shear modulus
    """
    xy, xz, zy = _get_pressures(energy_file)
    pressures = [np.mean([i,j,k]) for i,j,k in zip(xy,xz,zy)]
    volume = float(np.mean(trj.unitcell_volumes))
    volume *= 1e-27 * u.m**3
    temperature *= u.Kelvin

    GPa = 1e9*u.Pa
    u.define_unit("GPa", GPa)

    shear_bar, shear_std = _calc_mult(temperature, volume, pressures)
    shear_bar = shear_bar.in_units(u.GPa)
    shear_std = shear_std.in_units(u.GPa)

    return shear_bar, shear_std


def _calc_mult(temperature, volume, pressures):
    kb = 1.38E-23 * u.J / u.Kelvin
    constant = volume / (kb * temperature)
    p_list = list()
    p_unit = u.Pa
    for start_frame in np.linspace(0, len(pressures), num=5000, dtype=np.int):
        end_frame = start_frame + 2000
        if end_frame < len(pressures):
            p_chunk = pressures[start_frame:end_frame]
            print('\t\t\t...pressure {} to {}'.format(start_frame, end_frame))
            try:
                ensemble = np.mean([p**2 for p in p_chunk])
                ensemble *= p_unit ** 2
                total = constant * ensemble
                p_list.append(total)
            except TypeError:
                import pdb
                pdb.set_trace()
        else:
            continue
    p_bar = np.mean(p_list) * p_unit
    p_std = np.std(p_list) * p_unit
    return p_bar, p_std


def _get_pressures(energy_file):
    data = read_xvg(energy_file)
    data *= 0.986923 * u.atm
    data = data.to(u.Pa)
    xy = data[:,2]
    xz = data[:,3]
    zy = data[:,6]
    
    return xy, xz, zy
