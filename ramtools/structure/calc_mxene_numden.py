import numpy as np
import mdtraj as md
import warnings
import mdtraj.core.element as Element
import matplotlib.pyplot as plt
import os

def calc_number_density(coord_file, trj_file, bin_width, area, dim, box_range, data_path, resnames):
    """
    Calculate a 1-dimensional number density profile for each residue specifically for mxene water adsorption study

    Parameters
    ----------
    gro_file: str
        GROMACS '.gro' file to load 
    trj_file: str
        Trajectory to load
    top_file: str
        GROMACS '.top' file to load
    bin_width: int
        Width (nm) of numpy histogram bins
    dim: int
        Dimension to calculate number density profile (0,1 or 2)
    box_range: array
        Range of coordinates in 'dim' to evaluate
    data_path: str
        Path to save txt file out to
    resnames: dict
        Contains residues to look at

    
    Attributes
    ----------
    """

    first_frame = md.load_frame(trj_file, top=coord_file, index=0)

    open('{0}/resnames.txt'.format(data_path), 'w').close()

    for resname, restype in resnames.items():
        traj = md.load(trj_file, top=coord_file,
            atom_indices=first_frame.topology.select('name {}'
                .format(restype)))

        indices = [[at.index for at in compound.atoms]
            for compound in list(traj.topology.residues)]

        if 0 in [x.mass for x in
            [atom.element for atom in traj.topology.atoms]]:
            warnings.warn("mdtraj found zero mass, setting element to hydrogen", UserWarning)
            for atom in traj.topology.atoms:
                if atom.element in [Element.virtual, Element.virtual_site]:
                    atom.element = Element.hydrogen

        n_bins = 1 + round((box_range[1]-box_range[0]) / bin_width)
        x = np.histogram(traj.xyz[:, 1:, dim].reshape((-1, 1)),
            bins=np.linspace(box_range[0], box_range[1],
            num=n_bins))

        rho = x[0]
        bins = x[1]

        new_bins = (bins[:-1] + bins[1:]) / 2

        np.savetxt('{0}/{1}-number-density.txt'.format(data_path, resname),
            np.vstack([x[1][:-1]+np.mean(x[1][:2])-box_range[0],
            x[0]/(area*bin_width*(len(traj)-1))]).transpose())

        with open('{0}/resnames.txt'.format(data_path), "a") as myfile:
            myfile.write(resname + '\n')

def plot_mxene_numden(resnames, ylim, path, filename='number_density.pdf', shift=None):
    """
    function to plot number density profiles from txt files

    Paramters
    ---------
    resnames: list
        list of resnames to find corresponding txt files
    filename: str
        Name of PDF file to return
    shift: float, default=None
        Shift coordinates by this value
    
    Returns
    -------
    Number density profile in PDF format
    """
    def get_color(atom_name):
        """ Get color for matplotlib plot for each atom """
        color_dict = {
                'k': 'C0',
                'li': 'C1',
                'water': 'C2',
                'water_o': 'C9',
                'water_h': 'black',
                'O': 'C3',
                'OH': 'C4',
                'F': 'C7',
                'tam_N': 'black',
                'endC': 'red',
                'midC': 'blue',
                'branchC': 'green',
                }

        return color_dict[atom_name]

    def get_alpha(atom_name):
        """ Get color for matplotlib plot for each atom """
        alpha_dict = {
                'k': 1,
                'li': 1,
                'water': 1,
                'water_o': 1,
                'water_h': 1,
                'O': 0.2,
                'OH': 0.2,
                'F': 0.2,
                'tam_N': 1,
                'endC': 1,
                'midC': 1,
                'branchC': 1,
                }

        return alpha_dict[atom_name]

    fig, ax = plt.subplots()
    for f in resnames:
        data = np.loadtxt('{}/{}-number-density.txt'.format(path, f))
        if shift:
            ax.plot(data[:,0]-shift, data[:,1],
                    label=f,
                    color=get_color(f),
                    alpha=get_alpha(f)
                    )
        else:
            ax.plot(data[:,0], data[:,1],
                    label=f,
                    color=get_color(f),
                    alpha=get_alpha(f)
                    )

    plt.xlabel('Position on Surface (nm)')
    plt.ylabel('Number Density (nm^-3)')
    plt.ylim(ylim)
    plt.legend()
    plt.savefig(filename)

def calc_gmx_number_density(coord_file, trj_file, bin_width, area, dim, box_range, data_path, resnames, shift=None, chunk=None, unitcell_lengths=None, bond_array=None):
    """
    Calculate a 1-dimensional number density profile for each residue specifically for mxene water adsorption study

    Parameters
    ----------
    gro_file: str
        GROMACS '.gro' file to load 
    trj_file: str
        Trajectory to load
    top_file: str
        GROMACS '.top' file to load
    bin_width: int
        Width (nm) of numpy histogram bins
    dim: int
        Dimension to calculate number density profile (0,1 or 2)
    box_range: array
        Range of coordinates in 'dim' to evaluate
    data_path: str
        Path to save txt file out to
    resnames: dict
        Contains residues to look at
    shift: float, default=None
        shift coordinates by this value
    chunk: int, default=None
        Chunk of trajectory to consider
    unitcell_lengths: iterable, default=None
        Unitcell lengths of trajectory
    bond_array: iterable, default=None
        Array of bond indices

    
    Returns
    -------
    bins: Bins of histogram
    histogram_dict: dictionary of resname (key) and histogram values (value)
    """

    first_frame = md.load_frame(trj_file, top=coord_file, index=0)

    if os.path.exists(data_path):
        pass
    else:
        os.mkdir(os.path.join(os.path.abspath(os.getcwd()), data_path))

    open('{0}/resnames.txt'.format(data_path), 'w').close()
    
    histogram_dict = dict()
    traj = md.load(trj_file, top=coord_file)
    if unitcell_lengths:
        traj.unitcell_lengths = np.tile(unitcell_lengths, (traj.n_frames, 1))
        unitcell_vectors = np.array([[[unitcell_lengths[0], 0, 0],
                                     [0, unitcell_lengths[1], 0],
                                     [0, 0, unitcell_lengths[2]]]],
                                     dtype=np.float32,
        )
        traj.unitcell_vectors = np.tile(unitcell_vectors, (traj.n_frames, 1, 1))

    # Make molecules whole
    traj.make_molecules_whole(inplace=True, sorted_bonds=bond_array)

    for resname, restype in resnames.items():
        sub_traj = traj.atom_slice(traj.topology.select(
                   restype))
  
        if chunk:
            sub_traj = sub_traj[chunk:]

        indices = [[at.index for at in compound.atoms]
            for compound in list(sub_traj.topology.residues)]

        if 0 in [x.mass for x in
            [atom.element for atom in sub_traj.topology.atoms]]:
            warnings.warn("mdtraj found zero mass, setting element to hydrogen", UserWarning)
            for atom in sub_traj.topology.atoms:
                if atom.element in [Element.virtual, Element.virtual_site]:
                    atom.element = Element.hydrogen

        hist, bins = np.histogram(sub_traj.xyz[:, 1:, dim].reshape((-1, 1)),
            bins=np.linspace(box_range[0], box_range[1],
            num=1+round((box_range[1]-box_range[0])/bin_width)))
        bins_center = (bins[:-1] + bins[1:]) / 2

        if shift:
            np.savetxt('{0}/{1}-number-density.txt'.format(data_path, resname),
                np.vstack([bins_center - shift,
                hist/(area*bin_width*(len(sub_traj)-1))]).transpose())
        else:
            np.savetxt('{0}/{1}-number-density.txt'.format(data_path, resname),
                np.vstack([bins_center,
                hist/(area*bin_width*(len(sub_traj)-1))]).transpose())

        with open('{0}/resnames.txt'.format(data_path), "a") as myfile:
            myfile.write(resname + '\n')

        histogram_dict[resname] = hist / (area*bin_width*(len(sub_traj)-1))

    return bins_center, histogram_dict
