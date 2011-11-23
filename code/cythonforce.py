import numpy as np
cimport numpy as np
from math import sqrt


ctypedef float64
def norm(np.ndarray[np.float64_t, ndim=1] vec):
	""" Second norm """
	return sqrt(np.dot(vec,vec.conj()))

ctypedef np.float64_ dtype_t
def disp(np.ndarray[np.float64_t, ndim=1] one,
		np.ndarray[np.float64_t, ndim=1] two,
		np.ndarray[np.float64_t, ndim=1] offset):
    """ Displacement from one to two """
    cdef np.ndarray[np.float64_t, ndim=1] direct = one - two
    around = one + offset - two
    if norm(direct) <= norm(around):
        return direct
    else:
        return around