

def calc_conductivity(N,V,D_cat,D_an,q=1,T=300):
    kT = T * 1.3806488e-23
    q = q * 1.60218e-19

    cond = N / V * q ** 2 * (D_cat[0] + D_an[0]) / kT

    return cond
