import numpy as np
import os
import sys
import mdtraj as md


def find_atoms(gro_file, trj_file, top_file):
    trj = md.load(trj_file, top=gro_file)
    mins = [12.75,0,0]
    maxs = [15.46,0,0]
    indices = []
    for frame in trj.xyz:
        for atom_pos,atom in zip(frame,trj.topology.atoms):
            if mins[0] <= atom_pos[0] <= maxs[0]:
                indices.append(atom.index)
    #indices.append([atom_pos for atom_pos in frame if mins[0] <= atom_pos[0] <= maxs[0]])
    indices = np.array(indices)
    import pdb; pdb.set_trace()

def calc_number_density(gro_file, trj_file, top_file, bin_width, area,
        dim, box_range):
    """
    Calculate a 1-dimensional number density profile for each residue

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
    
    Attributes
    ----------
    """
    trj = md.load(trj_file, top=gro_file)
    resnames = np.unique([x.name for x in
               trj.topology.residues])
    mins = [12.75,0,0]
    maxs = [15.46,0,0]
    for resname in resnames:

        test = trj.topology.select('resname {}'.format(resname))
        """for frame in trj.xyz:
            for atom_pos in frame:
                if mins[0] <= atom_pos[0] <= maxs[0]: #what we want
                    import pdb; pdb.set_trace()"""
                    
        indices = [[atom.index for atom in compound.atoms] for compound in list(trj.atom_slice(test).topology.residues)]
        """indices = [[atom.index for atom in compound.atoms] if [mins[0]<=
            atom_pos[0] <= maxs[0]]] n"""
        import pdb; pdb.set_trace()
        com = np.array([md.compute_center_of_mass(trj.atom_slice(index)) for
            index in indices])
        x = np.histogram(com[:, 1:, dim].reshape((-1,1)),
            bins=np.linspace(box_range[0], box_range[1],
            num=1+round((box_range[1]-box_range[0])/bin_width)))

        np.savetxt('{}-number-density.txt'.format(resname),
                  np.vstack([x[1][:-1]+np.mean(x[1][:2])-box_range[0],
                  x[0]/(area*bin_width*(len(trj)-1))]).transpose())

        with open('resnames.txt', "a") as myfile:
            myfile.write(resname + '\n')


gro_file = 'sample.gro'
trj_file = 'sample.trr'
top_file = 'init.top'
box_range = [.92, 2.119] #2.119
dim = 1
bin_width = 0.01
area = 2.973 * 2.702
find_atoms(gro_file, trj_file, top_file)
