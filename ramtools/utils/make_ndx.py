import mdtraj as md
import numpy as np
import textwrap

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

    
    

