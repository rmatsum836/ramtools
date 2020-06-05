import numpy as np
import os
import sys
import mdtraj as md
from mtools.gromacs.gromacs import make_comtrj
import matplotlib as mpl
import matplotlib.pyplot as plt
from copy import deepcopy


def calc_number_density(trj, area, dim, box_range, n_bins,
        frame_range=None,maxs=None, mins=None, engine='gromacs'):
    """
    Calculate a 1-dimensional number density profile for each residue

    Parameters
    ----------
    trj: mdtraj trajectory
        trajectory from either GROMACS or LAMMPS
    bin_width: int
        Width (nm) of numpy histogram bins
    dim: int
        Dimension to calculate number density profile (0,1 or 2)
    box_range: array
        Range of coordinates in 'dim' to evaluate
    frame_range: Python range() (optional)
        Range of frames to calculate number density function over
    maxs: array (optional)
        Maximum coordinate to evaluate 
    mins: array (optional)
        Minimum coordinate to evalute
    engine: str, default = 'gromacs', additional option='lammps'
        Simulation engine of trajectory, determines unit conversions
    
    Attributes
    ----------
    """
    if engine == 'lammps':
        trj.unitcell_lengths *= 10
        trj.xyz[:,:,0] *= trj.unitcell_lengths[0][0]
        trj.xyz[:,:,1] *= trj.unitcell_lengths[0][1]
        trj.xyz[:,:,2] *= trj.unitcell_lengths[0][2]
        trj.xyz *= 10
    com_trj = make_comtrj(trj)
    resnames = np.unique([x.name for x in
               com_trj.topology.residues])
    rho_list = list()
    #bin_width = (box_range[1] - box_range[0]) / n_bins
    
    for resname in resnames:
        sliced = com_trj.topology.select('resname {}'.format(resname))
        trj_slice = com_trj.atom_slice(sliced)
        if frame_range:
            trj_slice = trj_slice[frame_range]
        for frame in trj_slice:
            if maxs is None:
                indices = [[atom.index for atom in compound.atoms]
                          for compound in
                          list(frame.topology.residues)]
            else:
                indices = np.intersect1d(
                          np.intersect1d(np.where(frame.xyz[-1, :, 0]
                              > mins[0]),
                          np.where(frame.xyz[-1, :, 0] < maxs[0])),
                          np.intersect1d(np.where(frame.xyz[-1, :, 1] 
                              > box_range[0]),
                          np.where(frame.xyz[-1, :, 1] < box_range[1])))

            if frame_range:
                if frame.time == frame_range[0]:
                    x = np.histogram(frame.xyz[0,indices,dim].flatten(), 
                        bins=n_bins, range=(box_range[0], box_range[1]))
                    rho = x[0]
                    bins = x[1]
                else:
                    rho += np.histogram(frame.xyz[0, indices, dim].
                            flatten(),bins=n_bins, range=(box_range[0],
                                box_range[1]))[0]
            else:
                if frame.time == 0:
                    x = np.histogram(frame.xyz[0,indices,dim].flatten(), 
                        bins=n_bins, range=(box_range[0], box_range[1]))
                    rho = x[0]
                    bins = x[1]
                else:
                    rho += np.histogram(frame.xyz[0, indices, dim].
                            flatten(),bins=n_bins, range=(box_range[0],
                                box_range[1]))[0]
            """else:
                rho += np.histogram(frame.xyz[0, indices, dim].flatten(), 
                          bins=n_bins, range=(box_range[0],
                          box_range[1]))[0]"""
        rho = np.divide(rho, trj_slice.n_frames * area *
                2 / n_bins)
        """print('With specified bins, rho = {}'.format(rho))
        test = np.divide(rho, trj_slice.n_frames * len(indices) 
                * area * 8  / 1000)
        print('With 1000 bins, rho = {}'.format(test))"""
        rho_list.append(rho)

    bin_list = bins[:-1]
    
    return(rho_list, bin_list)


def calc_area(trj):
    """
    Parameters
    ----------
    trj: mdtraj trajectory
    
    """
    clone = deepcopy(trj)
    clone.unitcell_lengths *= 10
    clone.xyz[:,:,0] *= clone.unitcell_lengths[0][0]
    clone.xyz[:,:,1] *= clone.unitcell_lengths[0][1]
    clone.xyz[:,:,2] *= clone.unitcell_lengths[0][2]
    clone.xyz *= 10

    resnames = np.unique([x.name for x in
           clone.topology.residues])
    resnames = np.delete(resnames, np.where(resnames=='RES'))
    sliced = clone.topology.select('resname {} or resname {}'.format(
       resnames[0], resnames[1]))
    trj_slice = clone.atom_slice(sliced)
    maxs = np.array([np.max(clone.xyz[:,:,x]) for x in range(3)])
    mins = np.array([np.min(clone.xyz[:,:,x]) for x in range(3)])
    intercal_region = maxs - mins
    area = np.prod(intercal_region) / 1000

    return area
