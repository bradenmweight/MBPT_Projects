import numpy as np
import subprocess as sp
from matplotlib import pyplot as plt
from math import factorial
import itertools


def get_globals():
    # Simulation Box
    global L, NL
    L = 10.0 # a.u.
    NL = 1000 # a.u.

    # External Potential
    global L_Well, V_well, N_SP_States
    L_Well = 4 # a.u.
    V_well = -10 # a.u.
    N_SP_States = 2

    # Particle Properties
    global m, N_electrons
    m = 1 # a.u.
    N_electrons = 2

def get_RGrid():
    return np.linspace( 0, L, NL )

def get_init_sp_wavefunctions():
    wfn_sp_dicts = []
    for j in range( 1, N_electrons+1 ):
        if ( j % 2 == 0 ):
            wfn_sp_dicts.append( {"sp_state":1, "spin":1} )
        else:
            wfn_sp_dicts.append( {"sp_state":1, "spin":-1} )
    
    return wfn_sp_dicts # Start all electrons in ground state SP orbital

def get_init_mb_wavefunctions( wfn_sp_dicts ):
    
    labels = [ wfn_sp_dicts[j]["sp_state"] for j in range(len(wfn_sp_dicts)) ]
    spins  = [ wfn_sp_dicts[j]["spin"] for j in range(len(wfn_sp_dicts)) ]

    wfn_mb_dictionaries = []

    for count, subset in enumerate( itertools.permutations(labels) ):
        subset = [ j for j in subset ]
        #print(subset)

        spin_labels = np.array(subset).astype(int) * -np.array(spins).astype(int)
        wfn_mb_dictionaries.append( { "spin-labels":(((count+1)%2)*2-1)*spin_labels, "sign":((count+1)%2)*2-1 } )

    #for d in range( len(wfn_mb_dictionaries) ):
    #    print( wfn_mb_dictionaries[d] )

    return wfn_mb_dictionaries

def get_external_potential( RGrid ):
    V = np.zeros(( NL ))
    for xi in RGrid:
        if ( 0 < xi and xi < L_Well  ):
            V[xi] = V_well
    return V

def get_SP_state( n, RGrid ):
    return np.sqrt(2/L_Well) * np.sin( n * np.pi * RGrid / L_Well )

def get_SP_energy( n ):
    return n**2 / ( 8 * np.pi * m * L_Well)

def get_init_rhos( wfn_sp_labels, RGrid ) :
    particle_rhos = []
    for j in wfn_sp_labels:
        state_j_left = np.conjugate(get_SP_state( j, RGrid ))
        state_j_right = get_SP_state( j, RGrid )
        particle_rhos.append( np.outer( state_j_left, state_j_right ) )
    return particle_rhos

def get_Coulomb_element_single_particle( indices, RGrid ):
    n,m,k,l = indices
    dR = RGrid[1]-RGrid[0]
    psi_n = get_SP_state( n, RGrid )
    psi_m = get_SP_state( m, RGrid )
    psi_k = get_SP_state( k, RGrid )
    psi_l = get_SP_state( l, RGrid )

    R_DIFF = np.abs( np.subtract.outer( RGrid, RGrid ) )
    R_DIFF[ np.diag_indices(NL) ] = 1.0
    V_int = 1 / R_DIFF
    V_int[ np.diag_indices(NL) ] = 0.0

    V_nmkl = 0
    for r1 in range( NL ):
        V_nmkl += np.sum( np.conjugate(psi_n)[r1] * np.conjugate(psi_m)[:] * V_int[r1,:] * psi_k[r1] * psi_l[:] )
    V_nmkl *= (RGrid[1] - RGrid[0])**2

    return V_nmkl

def save_1P_DMs( particle_rhos ):
    DM_DIR = "1P_DMs/"
    sp.call(f"mkdir -p {DM_DIR}",shell=True)

    for j in range(N_electrons):
        plt.imshow( particle_rhos[j], origin='lower' )
        plt.savefig(f"{DM_DIR}/density_matrix_{j}.jpg",dpi=300)

def get_all_Coulomb_elements( RGrid ):
    V = np.zeros(( N_SP_States, N_SP_States, N_SP_States, N_SP_States ))
    for n in range( 1, N_SP_States+1 ):
        for m in range( 1, N_SP_States+1 ):
            for k in range( 1, N_SP_States+1 ):
                for l in range( 1, N_SP_States+1 ):
                    V[n-1,m-1,k-1,l-1] = get_Coulomb_element_single_particle( [n,m,k,l], RGrid )
                    #print ( "n,m,k,l", n,m,k,l,V[n-1,m-1,k-1,l-1]  )
    return V

def main():
    get_globals()
    RGrid = get_RGrid()

    V_nmkl = get_all_Coulomb_elements( RGrid )

    wfn_sp_dicts = get_init_sp_wavefunctions()
    wfn_mb_dicts = get_init_mb_wavefunctions(wfn_sp_dicts)


    #print( "V_nmkl", V_nmkl )

    #particle_rhos = get_init_rhos( wfn_sp_labels, RGrid )
    #save_1P_DMs( particle_rhos )



    #V_nmkl = get_Coulomb_element( (1,2,3,4), RGrid )
    #print( V_nmkl )



if ( __name__ == "__main__" ):
    main()