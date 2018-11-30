

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
