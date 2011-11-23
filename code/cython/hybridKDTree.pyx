""" The goal here is to use Cython to speed up the scipy KDTree class """

import numpy as np
cimport numpy as np

dtype = np.float64

ctypedef np.float64_t dtype_t



def minkowski_distance_p_zero(np.ndarray[dtype_t, ndim=1] y,
                                float p=2.0):
    cdef float out
    if p==np.inf:
        out =  np.amax(np.abs(y-x),axis=-1)
    elif p==1:
        out =  np.sum(np.abs(y-x),axis=-1)
    else:
        out =  np.sum(np.abs(y-x)**p,axis=-1)
 
    return out

def minkowski_distance_p(np.ndarray[dtype_t, ndim=2] x,
                        np.ndarray[dtype_t, ndim=1] y,
                        float p=2.0):
    """Compute the pth power of the L**p distance between x and y

    For efficiency, this function computes the L**p distance but does
    not extract the pth root. If p is 1 or infinity, this is equal to
    the actual L**p distance.
    """
    #x = np.asarray(x)
    #y = np.asarray(y)
    cdef float out
    if p==np.inf:
        out =  np.amax(np.abs(y-x),axis=-1)
    elif p==1:
        out =  np.sum(np.abs(y-x),axis=-1)
    else:
        out =  np.sum(np.abs(y-x)**p,axis=-1)
 
    return out

def minkowski_distance(np.ndarray x,
                        np.ndarray y,
                        float p=2):
    """Compute the L**p distance between x and y"""
    #x = np.asarray(x)
    #y = np.asarray(y)
    if p==np.inf or p==1:
        return minkowski_distance_p(x,y,p)
    else:
        return minkowski_distance_p(x,y,p)**(1./p)


""" 
def minkowski_distance_p(x,y,p=2):
    #x = np.asarray(x)
    #y = np.asarray(y)
    if p==np.inf:
        return np.amax(np.abs(y-x),axis=-1)
    elif p==1:
        return np.sum(np.abs(y-x),axis=-1)
    else:
        return np.sum(np.abs(y-x)**p,axis=-1)
        
def minkowski_distance(x,y,p=2):
    x = np.asarray(x)
    y = np.asarray(y)
    if p==np.inf or p==1:
        return minkowski_distance_p(x,y,p)
    else:
        return minkowski_distance_p(x,y,p)**(1./p)

"""