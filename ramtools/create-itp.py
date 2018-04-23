def create_itp(itp_file, atoms, fx_const):
    f = open(itp_file, 'w')
    f.write('; position restraints built from "create_itp" function\n\n')
    f.write('[ position_restraints ]\n')
    f.write(';  i funct       fcx         fcy          fcz\n')
    count = range(1,atoms+1)
    for i in count:
        f.write('   {0}    1      {1}        {1}        {1}\n'.format(i,
            fx_const))
    f.close()

