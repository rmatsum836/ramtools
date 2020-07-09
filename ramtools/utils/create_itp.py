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
