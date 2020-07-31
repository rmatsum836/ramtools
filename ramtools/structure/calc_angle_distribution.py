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

def calc_water_angle(trj_file, gro_file, cutoff, system, dim=2):
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
    gro_file : GROMACS gro file
        MD coordinates to load
    cutoff : float
        Cutoff to analyze molecules in z-direction
    dim : int
        Dimension of surface vector
    """
    if dim == 0:
        normal_vector = [1, 0, 0]
    elif dim == 1:
        normal_vector = [0, 1, 0]
    else:
        normal_vector = [0, 0, 1]

    universe = mda.Universe(gro_file, trj_file)

    water_groups = universe.select_atoms('resname SOL')
    print("Unwrapping water molecules")
    transform = transformations.unwrap(water_groups)
    universe.trajectory.add_transformations(transform)
    print("Finished unwrapping water molecules")
    coordinates = [water_groups.positions for ts in universe.trajectory[5000:]]
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
                    np.array(normal_vector))

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
    plt.savefig(f'inverted_{system}_angles_{cutoff}.pdf')
    fig, ax = plt.subplots()
    #plt.bar(new_x, new_x_hist)
    plt.bar(new_x, y)
    plt.xlim((0, 181))
    plt.ylabel('Count')
    plt.xlabel('Angle (Deg)')
    plt.savefig(f'water_angles.pdf')
