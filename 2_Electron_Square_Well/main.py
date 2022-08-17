import numpy as np

def get_globals():
    # Simulation Box
    global L, NL
    L = 10.0 # a.u.
    NL = 1000 # a.u.

    # External Potential
    global L_Well, V_well, m
    L_Well = 4 # a.u.
    V_well = -1 # a.u.
    m = 1 # a.u.

def get_RGrid():
    return np.linspace( -L/2, L/2, NL )

def get_external_potential( RGrid ):
    V = np.zeros(( NL ))
    for xi in RGrid:
        if ( -L_Well/2 < xi and xi < L_Well/2  ):
            V[xi] = V_well
    return V

def get_SP_state( n, RGrid ):
    return np.sqrt(2/L) * np.sin( n * np.pi * RGrid / L )

def get_SP_energy( n ):
    return n**2 / ( 8 * np.pi * m * L_Well)

def get_Coulomb_element( indices, RGrid ):
    n,m,k,l = indices
    dR = RGrid[1]-RGrid[0]
    psi_n = get_SP_state( n, RGrid )
    psi_m = get_SP_state( m, RGrid )
    psi_k = get_SP_state( k, RGrid )
    psi_l = get_SP_state( l, RGrid )

    R_DIFF = np.subtract.outer( RGrid, RGrid )
    R_DIFF[ np.diag_indices(NL) ] = 1.0
    V_int = 1 / R_DIFF
    V_int[ np.diag_indices(NL) ] = 0.0


    # CHECK THIS ??? ~ BMW
    V_nmkl = 0
    for r1 in range( NL ):
        V_nmkl += np.sum( np.conjugate(psi_n)[r1] * np.conjugate(psi_m)[:] * V_int[r1,:] * psi_k[r1] * psi_l[:] )

    #V_nmkl = np.sum( np.conjugate(psi_n) * np.conjugate(psi_m) * \
    #            psi_k * psi_l  ) * dR**2
    return V_nmkl



def main():
    get_globals()
    RGrid = get_RGrid()
    V_nmkl = get_Coulomb_element( (1,2,3,4), RGrid )
    print( V_nmkl )



if ( __name__ == "__main__" ):
    main()