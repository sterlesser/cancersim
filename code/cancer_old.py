#standard imports, renamed for convenience
import scipy as sp
import pylab as py

#Initialization, let's create a hexagonal lattice
#so we will use a retangular unit cell, of width 1
#and height sqrt(3).  Each rectangle will have two guys

#make some epidermal cells
N = 10

XSIZE = 15
YSIZE = 5

#setup the top cells
epicells = []
for i in sp.arange(0,XSIZE,1):
    for j in sp.arange(0,YSIZE,1):
        if j:
            epicells.append([i,j*sp.sqrt(3)])
        epicells.append([i+0.5, (j+0.5)*sp.sqrt(3)])

epi_array = sp.array(epicells)

#setup the bottom cells
dermalcells = []
scale = 1.5
for i in sp.arange(0,XSIZE,step=scale):
    for j in sp.arange(0,YSIZE,step=scale):
        if j:
            dermalcells.append([i,-j*sp.sqrt(3)])
        dermalcells.append([(i+0.5*scale),-(j+0.5*scale)*sp.sqrt(3)])

derm_array = sp.array(dermalcells)

#setup the middle cells
middle_cells = []
for i in sp.arange(0,XSIZE,1):
    middle_cells.append([i,0])

mid_array = sp.array(middle_cells)

#put them all together into one array
cell_arr = sp.vstack( (epi_array,mid_array,derm_array ))
#define the cell types, 1:top, 2:middle, 3:bottom, 0:cancer
cell_types = [1]*len(epi_array) + [2]*len(mid_array) + [3]*len(derm_array)

#jiggle the cells
cell_arr += 0.1*sp.randn(*cell_arr.shape)


#use the Delaunay algorithm to set the initial links
from scipy.spatial import Delaunay
tri = Delaunay(cell_arr)
from collections import defaultdict

#store the links as a dictionary of sets
links = defaultdict(set)
for i,j,k in tri.vertices:
    links[i].update([j,k])
    links[j].update([i,k])
    links[k].update([i,j])


#compute the initial link lengths
from scipy.spatial.distance import euclidean
link_lengths = {}
for cell in links:
    for neigh in links[cell]:
        if cell > neigh:
            link_lengths[(neigh,cell)] = euclidean(cell_arr[cell],
                                                    cell_arr[neigh])
#Note that the link_lengths are in numerical order by index


#strength of the spring constants
spring_constants = [0.5,0.5,1.0,0.1]


def forces(cell_arr):
    """ Compute the force vector from the cell vectors """
    #cell_arr.shape = (cell_arr.size/2,2)
    forces = sp.zeros_like(cell_arr)
    #loop over all cells
    for cell,pos in enumerate(cell_arr):
        for neigh in links[cell]:
            if cell > neigh:
                #average the spring constants
                k = 0.5 * (spring_constants[cell_types[cell]]
                             + spring_constants[cell_types[neigh]])
                dist = euclidean(pos,cell_arr[neigh])
                diff = cell_arr[neigh] - pos
                L = link_lengths[(neigh,cell)]
                force = k*(dist-L)*diff/dist
                forces[cell] += force
                forces[neigh] -= force
 
    return forces


def energy(cell_arr):
    """ Compute the energy of the cell arrangement """
    #cell_arr.shape = (cell_arr.size/2,2)
    energy = 0.
    #loop over all cells
    for cell,pos in enumerate(cell_arr):
        for neigh in links[cell]:
            if cell > neigh:
                #average the spring constants
                k = 0.5 * (spring_constants[cell_types[cell]]
                             + spring_constants[cell_types[neigh]])
                
                dist = euclidean(pos,cell_arr[neigh])
                L = link_lengths[(neigh,cell)]
                energy += 0.5 * k * (dist-L)**2

    return energy


def link_energy(cell_arr,one,two):
    k = 0.5 * (spring_constants[cell_types[one]] 
                    + spring_constants[cell_types[two]])
    
    link_key = tuple(sorted((one,two)))
    dist = euclidean( cell_arr[one],cell_arr[two])
    L = link_lengths[link_key]

    return 0.5*k*(dist-L)**2





def add_bond(one,two,links):
    links[one].add(two)
    links[two].add(one)

def remove_bond(one,two,links):
    links[one].remove(two)
    links[two].remove(one)




#from scipy.optimize import fmin_bfgs, fmin_l_bfgs_b


test = cell_arr + 0.1*sp.randn(*cell_arr.shape)


from fire import FIRE
ans = FIRE(test,forces)






def plot_cells(cell_arr,fignum=1,*args,**kwargs):
    """ Plot the cells and the links """
    cell_colors = ['r','b','g','c']
    py.figure(fignum)
    py.clf()
    for cell in links:
        for neigh in links[cell]:
            if neigh>cell:
                data = cell_arr[[cell,neigh]]
                py.plot(data[:,0],data[:,1],
                    c=py.cm.jet( min(link_energy(cell_arr,cell,neigh)*30,1.) ),
                    alpha=0.8,zorder=1) 
                       
    py.scatter(cell_arr[:,0],cell_arr[:,1],
                    c=[cell_colors[i] for i in cell_types],
                    s=50,
                    zorder=10,
                    *args,**kwargs)
    py.axis('equal')