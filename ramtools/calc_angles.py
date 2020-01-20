import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis


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
    if theta > np.pi / 2:
        theta -= np.pi
    return theta

universe = MDAnalysis.Universe('nvt.tpr', 'nvt.xtc')
universe.select_atoms('type opls_774')

carbonyl_oxygen = universe.select_atoms('type opls_771')
ester_oxygen = universe.select_atoms('type opls_773')
methyl_carbon = universe.select_atoms('type opls_776')
carbonyl_carbon = universe.select_atoms('type opls_772')
ring_carbon = universe.select_atoms('type opls_774 opls_775')

atom_groups = {
    'Ring oxygen' : ester_oxygen,
    'Ring carbon' : ring_carbon,
    'Carbonyl oxygen' : carbonyl_oxygen,
    'carbonyl carbon' : carbonyl_carbon,
    'Methyl carbon' : methyl_carbon,
}

box = universe.select_atoms('all').bbox()

for name, atom_group in atom_groups.items():
    hist = np.zeros(400)
    min = 40
    max = 40
    for frame in universe.trajectory[::1]:
        tmp, bins = np.histogram(atom_group.positions[:, 2], bins=np.linspace(16.46, 119.15, len(hist) + 1))
        hist += tmp
        if np.min(bins) < min:
            min = np.min(bins)
        if np.max(bins) > max:
            max = np.max(bins)
    norm = box[1, 0] * box[1, 1] * universe.trajectory.n_frames * (bins[1] - bins[0])# * atom_group.n_atoms
    print(name, min, max)
    plt.plot(bins[1:] - np.mean(bins[:2]), hist / norm, label=name)

res_groups = [universe.select_atoms('resid {}'.format(r.resid)) for r in universe.residues if r.resname == 'PC']

running_list = list()

hist = np.zeros(400)
for frame in universe.trajectory[::100]:
    for i, atom_group in enumerate(res_groups):
        running_list.append(atom_group.center_of_mass()[2])
hist, bins = np.histogram(running_list, bins=np.linspace(16.46, 119.15, len(hist) + 1))

norm = box[1, 0] * box[1, 1] * universe.trajectory.n_frames * (bins[1] - bins[0])# * atom_group.n_atoms
plt.plot(bins[1:] - np.mean(bins[:2]), hist / norm, 'k--', label='COM')
print(np.max(hist))
plt.xlim((0, 20))
plt.ylim(ymin=0)
plt.xlabel(r'Position in channel, $\AA$')
plt.ylabel(r'Atomic density, $\AA^{-3}$')
plt.legend(loc=0)
plt.savefig('tmp.pdf')

ring_groups = [universe.select_atoms('type opls_773 opls_774 opls_775 and resid {}'.format(r.resid)) for r in universe.residues if r.resname == 'PC']

angles = list()
j =0
for ring_group in ring_groups:
    j += 1
    print(j)
    for frame in universe.trajectory[::1]:
        if ring_group.center_of_mass()[2] > 22:
            continue
        xyz = ring_group.positions
        tmp_A = []
        tmp_b = []
        for i in range(len(xyz[:, 0])):
            tmp_A.append([xyz[i, 0], xyz[i, 1], 1])
            tmp_b.append(xyz[i, 2])
        b = np.matrix(tmp_b).T
        A = np.matrix(tmp_A)
        fit = (A.T * A).I * A.T * b
        fit = np.array(fit.T)[0]
        angles.append(angle_between([fit[0], fit[1], -1], [0, 0, 1]) * 180 / np.pi)

plt.hist(angles, bins=180)
plt.xlim((-90, 90))
plt.savefig('angles.pdf')
