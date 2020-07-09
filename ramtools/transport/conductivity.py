import unyt as u

def calc_conductivity(N, V, D_cat, D_an, q=1, T=300):
    """ Calculate Nernst-Einstein Conductivity

    Parameters
    ----------
    N : int
        Number of ions
    V : float
        Volume of simulation box
    D_cat : float
        Diffusivity of cation in m^2/s
    D_an : float
        Diffusivity of anion in m^2/s
    q : float, default=1
        Charge of ions in element charge units
    T : float, default=300
        Temperature of system in Kelvin

    Returns
    -------
    cond : unyt.array
        Nernst-Einstein conductivity

    """
    D_cat *= u.m**2 / u.s
    D_an *= u.m**2 / u.s
    kT = T * 1.3806488e-23 * u.joule
    q *= u.elementary_charge
    q = q.to('Coulomb')

    cond = N / (V*kT) * q ** 2 * (D_cat + D_an)

    return cond
