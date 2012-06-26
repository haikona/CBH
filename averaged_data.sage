import time
import numpy as np

def averaged_data(L, filename, return_data=False):
    """
    INPUT:
    -- L: list of tuples (H, coeffs, data), where
       H: height
       coeffs: list of coefficients of elliptic curve giving that height
       data: the type of data that's being averaged, e.g., 2-Selmer, 3-Selmer
    -- filename: name of output file, passed in as a string

    OUTPUT:
    -- writes the averaged data to file
    -- (optional) returns numpy array of dimension len(list) x 2, 
       where the first column has heights, 
       and the second the averaged data up to that height
    """
    X = np.array([C[0] for C in L])
    V = np.array([C[2] for C in L])
    N = np.arange(1,X.shape[0]+1,dtype=np.float64)
    Y = np.cumsum(V)/N
    I = X[:-1] != X[1:]
    I = np.append(I,True)
    Z = np.vstack([X[I],Y[I]]).T
    np.savetxt(filename, Z)
    if return_data:
       return Z
