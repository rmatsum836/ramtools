import mdtraj as md
import numpy as np
import textwrap

def create_itp(itp_file, n_atoms, fx_const):
    """Create a GROMACS itp file for adding position restraints

    Parameters
    ----------
    itp_file : str
        Name of itp file to write to
    n_atoms : int
        Number of atoms in system
    fx_const : float
        Force constant for each atom in x, y, and z position
    """
    f = open(itp_file, 'w')
    f.write('; position restraints built from "create_itp" function\n\n')
    f.write('[ position_restraints ]\n')
    f.write(';  i funct       fcx         fcy          fcz\n')
    count = range(1,n_atoms+1)
    for i in count:
        f.write('   {0}    1      {1}        {1}        {1}\n'.format(i,
            fx_const))
    f.close()

def find_index_distance(gro_file, index_file, dist):
    """
    Find atoms in GROMACS based on distance criteria and append to an
    existing index file

    Parameters
    ----------
    gro_file: str
        gromacs .gro file used to look at atoms in system
    index_file: str
        gromacs index file to append to
    dist: float
        distance criteria for atoms to find
    """

    index_list = list()
    with open(gro_file, 'r') as f:
        for i,line in enumerate(f):
            if i < 3 or i > 9318:
                continue
            else:
                parse = line.split(' ')
                zdist = float(parse[-1].split('\n')[0])
                if zdist == dist:
                    index_list.append(int(parse[11]))

    with open(index_file, 'a+') as f:
        for i,index in enumerate(index_list):
            if i % 15 == 0 and i != 0:
                f.write('{} \n'.format(str(index)))
            else:
                f.write('{} '.format(str(index)))

def make_comtrj(trj):
    """Takes a trj and returns a trj with COM positions as atoms
    
    Note: Taken from mtools: https://github.com/mattwthompson/mtools/blob/master/mtools/gromacs/gromacs.py
    """

    comtop = md.Topology()
    coords = np.ndarray(shape=(trj.n_frames, trj.n_residues, 3))

    for j, res in enumerate(trj.topology.residues):
        comtop.add_atom(res.name, virtual_site, comtop.add_residue(res.name, comtop.add_chain()))
        res_frame = trj.atom_slice([at.index for at in res.atoms])
        coords[:, j, :] = md.compute_center_of_mass(res_frame)

    comtrj = md.Trajectory(xyz=coords,
                           topology=comtop,
                           time=trj.time,
                           unitcell_angles=trj.unitcell_angles,
                           unitcell_lengths=trj.unitcell_lengths)

    return comtrj


def make_ndx(top, selection1, selection2, output):
    """A substitute for gmx make_ndx that allow for more automation
    Inputs
    ______
    top: GROMACS .gro file 
    selection1: first resname 
    selection2: second resname
    output: Name of .ndx file

    Returns
    _______
    index: GROMACS .ndx file
        A .ndx file cointaining indices of a trajectory
    """

    topology = md.load(top).topology
    table, bonds = topology.to_dataframe()
    index_1  = table.loc[(table['resName']==selection1),
            'serial'].tolist()
    index_2 = table.loc[(table['resName']==selection2),
            'serial'].tolist()
    index_comb = np.append(index_1,index_2)
    out_file = open(output, "w")
    out_file.write('[ {0}-{1} ]\n'.format(selection1, selection2))
    atoms = '{}\n'.format('    '.join(str(value) for value in index_comb))
    out_file.write(textwrap.fill(atoms,80))
    out_file.write('\n')
    out_file.close()
    return(out_file)

    
