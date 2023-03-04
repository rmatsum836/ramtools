from mtools.post_process import calc_msd
import mdtraj as md
import json

def _run_overall(trj, mol):
     D, MSD, x_fit, y_fit = calc_msd(trj)
     return D, MSD

def _run_multiple(trj, mol):
    D_pop = list()
    for start_frame in np.linspace(0, 1001, num=200, dtype=np.int):
        end_frame = start_frame + 200
        if end_frame < 1001:
            chunk = trj[start_frame:end_frame]
            print('\t\t\t...frame {} to {}'.format(start_frame, end_frame))
            try:
                D_pop.append(calc_msd(chunk)[0])
            except TypeError:
                import pdb
                pdb.set_trace()
        else:
            continue
    D_bar = np.mean(D_pop)
    D_std = np.std(D_pop)
    return D_bar, D_std
