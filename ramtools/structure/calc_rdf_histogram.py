import numpy as np
import mdtraj as md
import itertools
import matplotlib.pyplot as plt


def calc_rdf_histo(trj, pairs, r_range=None, bin_width=0.005, n_bins=None,
        periodic=True, opt=True):
    """
    This is taken straight from mdtraj just to get the
    histogram of the RDF
    """
    if r_range is None:
        r_range = np.array([0.0, 1.0])
    r_range = md.utils.ensure_type(r_range, dtype=np.float64, ndim=1, name='r_range',
                          shape=(2,), warn_on_cast=False)
    if n_bins is not None:
        n_bins = int(n_bins)
        if n_bins <= 0:
            raise ValueError('`n_bins` must be a positive integer')
    else:
        n_bins = int((r_range[1] - r_range[0]) / bin_width)

    distances = md.geometry.distance.compute_distances(trj, pairs, periodic=periodic, opt=opt)
    g_r, edges = np.histogram(distances, range=r_range, bins=n_bins)

    return g_r, edges

"""
GUIDE TO MXENE SYSTEM ATOM TYPES

Example: trj.topology.select('name "9"')

1 6.9412 Lithium
2 47.8671 Titanium
3 47.8671 Titanium
4 12.0108 MXene Carbon
5 15.9994 MXene OH oxygen
6 1.00795 MXene H hydrogen
7 18.9984 MXene Fluorine
8 15.9994 MXene Oxygen
9 15.9994 Water Oxygen
10 1.00795 Water Hydrogen

"""

trj = md.load('md_step.dcd', top='data.gro')

# Get combinations of atom names to get pairs
types = {'li': '"1"', 'water': '"9" "10"', 'oh': '"5" "6"', 'f': '"7"', 'o': '"8"'}
#combos = itertools.combinations(['"1"', '"9" "10"', '"5" "6"', '"7"', '"8"'], 2)
combos = itertools.combinations([i for i in types.values()], 2)
names = itertools.combinations([i for i in types], 2)
for combo,name in zip(combos, names):
    fig, ax = plt.subplots()
    pairs = trj.topology.select_pairs(selection1='name {}'.format(combo[0]),
            selection2='name {}'.format(combo[1]))
    g_r, edges = calc_rdf_histo(trj, pairs=pairs)
    width = 0.7 * (edges[1] - edges[0])
    center = (edges[:-1] + edges[1:]) / 2
    ax.bar(center, g_r, width=width)
    plt.xlabel('distance (nm)')
    plt.ylabel('bins')
    plt.savefig('rdf_hist_{}_{}.pdf'.format(name[0], name[1]))
