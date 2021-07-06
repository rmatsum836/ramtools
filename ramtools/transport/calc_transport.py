import numpy as np
import mdtraj as md
import unyt as u

from ramtools.utils.read_files import read_xvg
from scipy import stats


def calc_ne_conductivity(N, V, D_cat, D_cat_error, D_an, D_an_error, q=1, T=300):
    """ Calculate Nernst-Einstein Conductivity

    Parameters
    ----------
    N : int
        Number of ions
    V : float
        Volume of simulation box in nm^3
    D_cat : float
        Diffusivity of cation in m^2/s
    D_cat_error : float
        The error of diffusivity of cation in m^2/s
    D_an : float
        Diffusivity of anion in m^2/s
    D_an_error : float
        The error of diffusivity of anion in m^2/s
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
    cond_error = N / (V*kT) * q ** 2 * np.sqrt(np.square(D_cat_error) + np.square(D_an_error))
    return cond, cond_error

def calc_eh_conductivity(trj_file, gro_file, V, cat_resname, an_resname, chunk=200, cut=0, q=1, T=300,
                         skip=100):
    """ Calculate Einstein-Helfand conductivity
    Parameters
    ----------
    trj_file : GROMACS trr or xtc file
        GROMACS trajectory
    gro_file : GROMACS gro file
        GROMACS coordinate file
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
    overall_avg = list()
    trj = md.load(trj_file, top=gro_file)

    n_frames = trj.n_frames
    for outer_chunk in range(2000, n_frames, 2000):
        running_avg = np.zeros(chunk)
        trj_outer_chunk = trj[outer_chunk-2000:outer_chunk]
        for i, start_frame in enumerate(np.linspace(0, 2000, num=500, dtype=np.int)):
            end_frame = start_frame + chunk
            if end_frame < 2000:
                trj_chunk = trj_outer_chunk[start_frame:end_frame]
                if i == 0:
                    trj_time = trj_chunk.time
                trj_slice = trj_chunk.atom_slice(trj_chunk.top.select(f'resname {cat_resname} {an_resname}'))
                total_n_il = trj_slice.n_atoms
                q_pos = np.array([q] * int(total_n_il / 2))
                q_neg = -1 * q_pos
                q_list = np.append(q_pos, q_neg)
                M = trj_slice.xyz.transpose(0, 2, 1).dot(q_list)
                running_avg += [np.linalg.norm((M[i] - M[0])) ** 2 for i in range(len(M))]
        y = running_avg / i
        overall_avg.append(y)

    x = (trj_time - trj_time[0]).reshape(-1)
    sigma_all = list()
    kB = 1.38e-23 * u.joule / u.Kelvin
    V *= u.nm ** 3
    T *= u.Kelvin
    for i in range(len(overall_avg)):
        slope, intercept, r_value, p_value, std_error = stats.linregress(x[cut:], overall_avg[i][cut:])
        sigma = slope * (u.elementary_charge * u.nm) ** 2 / u.picosecond / (6 * V * kB * T)
        seimens = (u.Ω) ** (-1)
        sigma = sigma.to(seimens / u.meter)

        sigma_all.append(sigma)
    std = np.std(sigma_all)
    sigma_new = np.average(sigma_all)
    y_new = np.average(overall_avg, axis =0)

    return x, y_new, sigma_new, std


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
