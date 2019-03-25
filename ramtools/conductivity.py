

def calc_conductivity(q=1,N,V,D_cat,D_an):
    kT = T * 1.3806488e-23
    q = q * 1.60218e-19

    cond = N / V * q ** 2 * (D_cat + D_an) / kT

    return cond
