import numpy as np
import mdtraj as md
import parmed as pmd
import warnings
import mdtraj.core.element as Element
import matplotlib.pyplot as plt

def calc_conductivity(traj, N_pairs, T, V, D_cat, D_an, q_cat,
        q_an):
    """
    Calculate the nernt-einstein conductivity of a bulk simulation

    Parameters
    ----------
    traj: mdtraj trajectory
    N_pairs: Number of ion pairs
    T: Temperature of system
    V: Volume of system
    D_cat: diffusivity of cation
    D_an: diffusivity of anion
    q_cat: charge of cation
    q_an: charge of anion

    """
    kb = 1.38E-23 # m^2kgs^-2K-1 Boltzmann constant
    sigma = (N_pair)/(V*kb*T)*((q_cat**2*D_cat)+(q_an**2*D_an))

