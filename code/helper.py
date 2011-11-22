import scipy as sp

from math import sqrt
from scipy import dot

from scipy.linalg import norm as builtinnorm

########################################################
### Helper functions ###################################
########################################################

def normold(vec):
    """ Normalize a vector """
    return sp.sqrt(sp.sum(vec**2))

def norm(vec):
	""" Second norm """
	return sqrt(dot(vec,vec.conj()))

def norm3(vec):
	""" Third norm """
	return builtinnorm(vec)

def unitize(vec):
    """ Construct a unit vector from a vector """
    return vec/norm(vec)

