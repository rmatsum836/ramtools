import unyt as u

def calc_conductivity(N,V,D_cat,D_an,q=1,T=300):
    D_cat *= u.m**2 / u.s
    D_an *= u.m**2 / u.s
    kT = T * 1.3806488e-23 * u.joule
    q *= u.elementary_charge
    q = q.to('Coulomb')

    cond = N / (V*kT) * q ** 2 * (D_cat + D_an)

    return cond
