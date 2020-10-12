import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis import transformations
from scipy.ndimage.measurements import center_of_mass
from scipy import stats

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

    >>> angle_between((1, 0, 0), (0, 1, 0))
    1.5707963267948966
    >>> angle_between((1, 0, 0), (1, 0, 0))
    0.0
    >>> angle_between((1, 0, 0), (-1, 0, 0))
    3.141592653589793
    """
    v1 = v1 / np.linalg.norm(v1)
    v2 = v2 / np.linalg.norm(v2)
    theta = np.arccos(np.clip(np.dot(v1, v2), -1.0, 1.0))

    return theta

def s_order_parameter(angles):
    """ Calculate angle order parameter
    
    Parameters
    ----------
    angles : list
        Angle between water and surface in radians
    """
    angle_average = np.mean(np.cos(angles)**2)
    s = (3*angle_average-1) / 2

    return s

def calc_water_angle(trj_file, gro_file, cutoff, dim=2, filepath=''):
    """ Calculate angle distribution between a water molecule vector and normal of
    a surface

    Water vector:

        ^
        |
        |

        O
       / \
      H   H

    Parameters
    ----------
    trj_file : trajectory file
        MD trajectory to load
    gro_file : Coordinate file
        MD coordinates to load.  MOL2 file is preferred as it contains bond information.
    cutoff : float
        Cutoff to analyze molecules in z-direction (angstroms)
    dim : int
        Dimension of surface vector
    """
    if dim == 0:
        normal_vector = [1, 0, 0]
    elif dim == 1:
        normal_vector = [0, 1, 0]
    else:
        normal_vector = [0, 0, 1]

    trj_str = f'{filepath}/{trj_file}'
    gro_str = f'{filepath}/{gro_file}'

    universe = mda.Universe(gro_str, trj_str)

    water_groups = universe.select_atoms('resname SOL')
    print("Unwrapping water molecules")
    transform = transformations.unwrap(water_groups)
    universe.trajectory.add_transformations(transform)
    print("Finished unwrapping water molecules")
    coordinates = [water_groups.positions for ts in universe.trajectory]
    angles = list()
    radians = list()
    print("Starting to analyze vectors ... ")
    for frame_num, frame in enumerate(coordinates):
        for idx in np.arange(3,len(frame)+3,3):
            xyz = frame[idx-3:idx]

            if xyz[0][dim] > cutoff:
                continue
            # Get midpoint of hydrogens
            fit = [(xyz[1][i]+xyz[2][i])/2 for i in range(3)]
            # Draw vector of oxygen going through hydrogen midpoint
            vector = [xyz[0][i] - fit[i] for i in range(3)]

            angle = angle_between(np.array([vector[0], vector[1], vector[2]]),
                    np.array(normal_vector)) * (180 / np.pi)

            angle_in_radians = angle * np.pi / 180
            radians.append(angle_in_radians)
            angles.append(angle)


    y, x = np.histogram(angles, bins=180, density=True, range=(0.0, 180.0))
    new_x = list()
    for idx in range(180):
        mid = idx + 0.5
        new_x.append(mid)
    new_x_hist = y / np.sin((np.array(new_x) * np.pi / 180))
    fig, ax = plt.subplots()
    plt.plot(new_x, y)
    plt.xlim((0, 181))
    plt.ylabel('Count')
    plt.xlabel('Angle (Deg)')
    fig, ax = plt.subplots()
    #plt.bar(new_x, new_x_hist)
    plt.bar(new_x, y)
    plt.xlim((0, 181))
    plt.ylabel('Count')
    plt.xlabel('Angle (Deg)')
    plt.savefig(f'{filepath}/water_angles.pdf')

def calc_water_order_parameter(trj_file, gro_file, cutoffs, bin_width=0.2, shift=True, dim=2, filepath=''):
    """ Calculate the order parameter between a water molecule vector and normal of
    a surface
    DOI: 10.1021/la0347354

    Water vector:

        ^
        |
        |

        O
       / \
      H   H

    Parameters
    ----------
    trj_file : trajectory file
        MD trajectory to load
    gro_file : Coordinate file
        MD coordinates to load.  MOL2 file is preferred as it contains bond information.
    cutoffs : list or tuple
        Dimensions of slitpore to consider (angstroms)
    bin_width : float, default = 0.2 angstroms
        bin width of histogram
    shift : boolean, default=True
        Shift center to 0 if True
    dim : int
        Dimension of surface vector
    filepath : str, default=''
        filepath to save plot to
    """
    if dim == 0:
        normal_vector = [1, 0, 0]
    elif dim == 1:
        normal_vector = [0, 1, 0]
    else:
        normal_vector = [0, 0, 1]

    trj_str = f'{filepath}/{trj_file}'
    gro_str = f'{filepath}/{gro_file}'

    universe = mda.Universe(gro_str, trj_str)

    water_groups = universe.select_atoms('resname SOL')
    print("Unwrapping water molecules")
    transform = transformations.unwrap(water_groups)
    universe.trajectory.add_transformations(transform)
    print("Finished unwrapping water molecules")
    coordinates = [water_groups.positions for ts in universe.trajectory]
    print("Starting to analyze vectors ... ")
    angle_position = dict()
    for frame_num, frame in enumerate(coordinates):
        for idx in np.arange(3,len(frame)+3,3):
            xyz = frame[idx-3:idx]

            if xyz[0][dim] < cutoffs[0] and xyz[0][dim] > cutoffs[1]:
                continue
            # Get midpoint of hydrogens
            fit = [(xyz[1][i]+xyz[2][i])/2 for i in range(3)]
            # Draw vector of oxygen going through hydrogen midpoint
            vector = [xyz[0][i] - fit[i] for i in range(3)]

            angle = angle_between(np.array([vector[0], vector[1], vector[2]]),
                    np.array(normal_vector)) * (180 / np.pi)

            angle_in_radians = angle * np.pi / 180
            angle_position[xyz[0][dim]] = angle_in_radians

    distance_dict = dict()
    for dis in np.arange(cutoffs[0], cutoffs[1], step=bin_width):
        distance_dict[(dis+(dis+bin_width))/2] = list()
        for pos, angle in angle_position.items():
            if pos > dis and pos < (dis+bin_width):
                distance_dict[(dis+(dis+bin_width))/2].append(angle)

    s_order_dict = dict()
    for dis, angles in distance_dict.items():
        s_order_dict[dis] = s_order_parameter(angles)

    if shift:
        new_bins = list(s_order_dict.keys())
        shift_value = (cutoffs[1] - cutoffs[0]) / 2 + cutoffs[0]
        new_bins = [(bi-shift_value) for bi in new_bins]

    fig, ax = plt.subplots()
    # Divide by 10 to go to nm
    if shift:
        final_bins = [i/10 for i in new_bins]
        plt.plot(final_bins, s_order_dict.values())
    else:
        final_bins = [i/10 for i in s_order_dict.keys()]
        plt.plot(final_bins, s_order_dict.values())
    plt.xlabel('Distance (nm)')
    plt.ylabel('S')
    if shift:
        plt.xlim((-1, 1))
    plt.ylim((-0.5, 0.25))
    plt.savefig(f'{filepath}/s_order.pdf')
    np.savetxt(f'{filepath}/s_order.txt',
               np.transpose(np.vstack([final_bins, list(s_order_dict.values())])))
