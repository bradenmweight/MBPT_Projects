import numpy as np

def get_globals():
    global L
    L = 1.0 # a.u.

def get_SP_state( n, RGrid ):
    return np.sqrt(2/L) * np.sin( n * np.pi * RGrid / L )
    