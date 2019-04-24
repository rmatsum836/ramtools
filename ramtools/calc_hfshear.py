import numpy as np
import mdtraj as md
import unyt as u


def read_xvg(fname):
    data=[]
    with open(fname) as f:
        for line in f:
            # Lines with metadata or comments start with #, @
            if not line.startswith(("@","#")):
                data.append(np.array([float(s) for s in line.split()]))
    data = np.vstack(data)
    return data


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


def calc_hfshear(energy_file, trj, temperature):
    xy, xz, zy = _get_pressures(energy_file)
    pressures = [np.mean([i,j,k]) for i,j,k in zip(xy,xz,zy)]
    volume = float(np.mean(trj.unitcell_volumes))
    volume *= 1e-27 * u.m**3
    temperature *= u.Kelvin

    shear_bar, shear_std = _calc_mult(temperature, volume, pressures)
    shear_bar *= 1e-9
    shear_std *= 1e-9

    return shear_bar, shear_std
