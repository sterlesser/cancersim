""" The goal here is to use Cython to speed up the scipy KDTree class """

import numpy as np
cimport numpy as np

import cython



from math import sqrt

dtype = np.float64

ctypedef np.float64_t dtype_t


@cython.boundscheck(False)
@cython.wraparound(False)
def disp_func(np.ndarray[dtype_t, ndim=1] x1,
            np.ndarray[dtype_t, ndim=1] x2,
            float XSIZE ):
    cdef np.ndarray[dtype_t, ndim=1] disp1, disp2, disp3
    disp1 = x1 - x2

    cdef float norm1, norm2, norm3
    norm1 = norm(disp1)

    if norm1 < 3.0:
        return disp1

    disp2 = x1 + XSIZE - x2
    norm2 = norm(disp2)

    disp3 = x1 - XSIZE - x2
    norm3 = norm(disp3)

    if norm1 <= norm2 and norm1 <= norm3:
        return disp1
    elif norm2 <= norm1 and norm2 <= norm3:
        return disp2
    else:
        return disp3


@cython.boundscheck(False)
@cython.wraparound(False)
def norm(np.ndarray[dtype_t, ndim=1] x):
    cdef float size
    size = 0.0
    for i in range(x.shape[0]):
        size += x[i]**2
    size = sqrt(size)
    return size

@cython.boundscheck(False)
@cython.wraparound(False)
def force_func_basal(np.ndarray[dtype_t, ndim=1] x1,
                     np.ndarray[dtype_t, ndim=1] x2,
                     float basalstrength,
                     float basalcutoff, 
                     float XSIZE ):
    #We have one basal cell
    cdef np.ndarray[dtype_t, ndim=1] disp
    disp = disp_func(x1,x2,XSIZE)

    cdef float mod_disp 
    mod_disp = norm(disp)

    cdef np.ndarray[dtype_t, ndim=1] force
    force = np.zeros(2, float)

    if mod_disp <= basalcutoff:
        force = 2 * basalstrength**4 * ( 2 * basalcutoff**2 - 3 * basalcutoff * mod_disp + mod_disp**2 )/( basalcutoff**2 * mod_disp**6 ) * disp

    return force






@cython.boundscheck(False)
@cython.wraparound(False)
def force_func_hertz(np.ndarray[dtype_t, ndim=1] x1,
                                np.ndarray[dtype_t, ndim=1] x2,
                                float r1,
                                float r2,
                                float a,
                                float XSIZE):
    """ Compute the Hertz force between two cells, takes
            x1 - pos of first cell
            x2 - pos of second cell
            r1 - radius of first cell
            r2 - radius of second cell
            a - strength of force """
    cdef np.ndarray[dtype_t, ndim=1] disp
    disp = disp_func(x1,x2,XSIZE)

    cdef float mod_disp 
    mod_disp = norm(disp)

    cdef float delta
    delta=(r1+r2)-mod_disp

    cdef np.ndarray[dtype_t, ndim=1] force
    force = np.zeros(2, float)

    if delta > 0:
        force = sqrt(r1*r2/(r1+r2)) * a * delta**1.5*disp/mod_disp

    return force