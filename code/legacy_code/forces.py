
from helper import norm, unitize
import scipy as sp


###################################
# Forces ##########################
###################################

def force_func(cell1,cell2,a,xi):
    """ the native force function between two positions """
    x1 = cell1.pos
    x2 = cell2.pos
    disp = x1 - x2
    mod_disp = norm(disp)
    force = 2 * a**4 * ( 2 * xi**2 - 3 * xi * mod_disp + mod_disp**2 )/( xi**2 * mod_disp**6 ) * disp
    
    return force

def force_func2(cell1,cell2,a,xi):
    """ the native force function between two positions """
    x1 = cell1.pos
    x2 = cell2.pos
    r1 = cell1.radius
    r2 = cell2.radius
    disp = x1 - x2
    mod_disp = norm(disp)
    a1=a*(r1+r2)
    xi1=xi*(r1+r2)
    force = 2 * a1**4 * ( 2 * xi1**2 - 3 * xi1 * mod_disp + mod_disp**2 )/( xi1**2 * mod_disp**6 ) * disp
    
    return force

def force_func_hertz(cell1,cell2,a,xi):
     """ the Hertz force between two cells """
     x1 = cell1.pos
     x2 = cell2.pos
     r1 = cell1.radius
     r2 = cell2.radius
     disp = x1 - x2
     mod_disp = norm(disp)
     delta=(r1+r2)-mod_disp
     if delta > 0.0:
         force = a*delta**1.5*disp/mod_disp
     else:
         force= 0.0

     return force
