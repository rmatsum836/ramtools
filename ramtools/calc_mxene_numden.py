import numpy as np
import mdtraj as md
import warnings
import mdtraj.core.element as Element

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

        x = np.histogram(traj.xyz[:, 1:, dim].reshape((-1, 1)),
            bins=np.linspace(box_range[0], box_range[1],
            num=1+round((box_range[1]-box_range[0])/bin_width)))

        np.savetxt('{0}/{1}-number-density.txt'.format(data_path, resname),
            np.vstack([x[1][:-1]+np.mean(x[1][:2])-box_range[0],
            x[0]/(area*bin_width*(len(traj)-1))]).transpose())

        with open('{0}/resnames.txt'.format(data_path), "a") as myfile:
            myfile.write(resname + '\n')
