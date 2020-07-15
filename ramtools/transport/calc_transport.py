import numpy as np
import mdtraj as md
import unyt as u

from ramtools.utils.read_files import read_xvg
from scipy import stats


def calc_ne_conductivity(N, V, D_cat, D_an, q=1, T=300):
    """ Calculate Nernst-Einstein Conductivity

    Parameters
    ----------
    N : int
        Number of ions
    V : float
        Volume of simulation box in nm^3
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
    V *= u.nm ** 3
    V = V.to('m**3')

    cond = N / (V*kT) * q ** 2 * (D_cat + D_an)

    return cond

def calc_eh_conductivity(trj_file, gro_file, N, V, cat_resname, an_resname, chunk=200, q=1, T=300,
        skip=100):
    """ Calculate Einstein-Helfand conductivity
    Parameters
    ----------
    trj_file : GROMACS trr or xtc file
        GROMACS trajectory
    gro_file : GROMACS gro file
        GROMACS coordinate file
    N : int
        Number of ions
    V : float
        Volume of simulation box
    cat_resname : str
        Residue name of cation
    an_resname : str
        Residue name of anion
    q : float, default=1
        Charge of ions in element charge units
    T : float, default=300
        Temperature of system in Kelvin
    skip : int, default=100
        Number of frames in trajectory to skip

    Returns
    -------
    cond : unyt.array
        Einstein-Helfand conductivity
    """

    running_avg = np.zeros(chunk)
    for i,trj in enumerate(md.iterload(trj_file, top=gro_file, chunk=chunk, skip=skip)):
        if i == 0:
            trj_time = trj.time
        if trj.n_frames != chunk:
            continue
        try:
            trj = trj.atom_slice(trj.top.select(f'resname {cat_resname} {an_resname}'))
        except:
            print("Not slicing trajectory")
        M = dipole_moments_md(trj, q)
        running_avg += [np.linalg.norm((M[i] - M[0]))**2 for i in range(len(M))]

    x = (trj_time - trj_time[0]).reshape(-1)
    y = running_avg / i

    # TODO: Find where slope becomes linear
    slope, intercept, r_value, p_value, std_error = stats.linregress(
            x, y)

    kB = 1.38e-23 * u.joule / u.Kelvin
    V *= u.nm ** 3
    T *= u.Kelvin

    sigma = slope * (u.elementary_charge * u.nm) ** 2 / u.picosecond / (6 * V * kB * T)
    seimens = u.second ** 3 * u.ampere ** 2 / (u.kilogram * u.meter ** 2)
    sigma = sigma.to(seimens / u.meter)

    return sigma


def calc_hfshear(energy_file, trj, temperature):
    """ Calculate High-Frequency shear modulus of an MDTraj trajectory

    Parameters
    ----------
    energy_file : str
        GROMACS .xvg file, created by running "gmx energy -f energy.edr -vis"
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

    shear_bar, shear_std = _calc_mult(temperature, volume, pressures)
    shear_bar = shear_bar.in_units(u.GPa)
    shear_std = shear_std.in_units(u.GPa)

    return shear_bar, shear_std

def dipole_moments_md(traj, charges):
    local_indices = np.array([(a.index, a.residue.atom(0).index) for a in traj.top.atoms], dtype='int32')
    local_displacements = md.compute_displacements(traj, local_indices, periodic=False)

    molecule_indices = np.array([(a.residue.atom(0).index, 0) for a in traj.top.atoms], dtype='int32')
    molecule_displacements = md.compute_displacements(traj, molecule_indices, periodic=False)

    xyz = local_displacements + molecule_displacements

    moments = xyz.transpose(0, 2, 1).dot(charges)

    return moments


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
