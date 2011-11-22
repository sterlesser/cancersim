#Cancer Sim

import scipy as sp
import pylab as py

from collections import defaultdict

from scipy.spatial.distance import euclidean
from math import pow
from scipy.spatial import Delaunay
from scipy.spatial import KDTree
import random
import time

XSIZE = 20
YSIZE = 20

########################################################
### Helper functions ###################################
########################################################

def norm(vec):
    """ Normalize a vector """
    return sp.sqrt(sp.sum(vec**2))

def unitize(vec):
    """ Construct a unit vector from a vector """
    return vec/norm(vec)


########################################################
### Cell Stuff  ########################################
########################################################

class CellType(object):
    """ A holder for some cell type information
        
        Initialization / Attributes:
            * name : string name of cell type 
            * type_ind : an int for the cell type  
            * k : spring constant (float) 
            * L : natural length of bonds (float) 
            * color='b' : set the color the plots
    """
        
    def __init__(self,name,type_ind,k,L,color='b'):
        self.name = name
        self.type_ind = type_ind 
        self.k = k 
        self.L = L  
        self.color = color
        
    def __repr__(self):
        return "<CellType: {0.name}, type_ind={0.type_ind}, k={0.k}, L={0.L}>".format(self)

class Cell(CellType):
    """ A cell object
        
        Initialization / Attributes:
            * pos: a position array
            * type: a type declaration, instance of CellType class
            * vel: a velocity arrays
            * index: an integer index for finding the cell in a big position array
    """
    def __init__(self,pos,cell_type,index=0):
        self.pos = sp.array(pos)
        self.type = cell_type
        self.vel = sp.zeros_like(self.pos)
        self.index = index
        self.radius = self.type.L/2.0

    def __repr__(self):
        return "<Cell type:{0.type.name} pos:{0.pos} ID:{1}>".format(self,id(self))


class GhostCell(Cell):
    """ A ghost cell, which is used to implement periodic boundary conditions.
        A ghost cell is linked to its parent.
        
        Initialization:
            * original: parent cell
            * sign: either 1 or -1
            * xsize: size of box, size of offset
    """
    def __init__(self,original,sign=1,xsize=XSIZE):
        self.offset = sign* sp.array([XSIZE,0.])
        self.sign = sign
        self.original = original
        self.type = original.type
        self.vel = original.vel

        self.pos = self.original.pos + self.offset
        
    def update(self):
        """ Update the position of the GhostCell based on the offset and sign """
        self.pos = self.original.pos + self.offset

    def __repr__(self):
        return "<GhostCell: pos:{0.pos} of cell {0.original}>".format(self)

########################################################
### Link Stuff  ########################################
########################################################


class Link(object):
    """ A link object keeps a link between two cells.
    
        Initialization:
            * one: cell one
            * two: cell two
            * L: resting length of link_disp
            * k: strength of link_disp
            * xsize: size of box (to handle periodic boundary conditions)

        Properties:
            * disp - displacement of spring
            * extension - compute the current link extension
            * energy - compute the current link energy
            * force - force of spring
    """
    def __init__(self,one,two,L=None,k=None,xsize=XSIZE):
        self.one = one
        self.two = two
        self.xsize=xsize
        self.offset = sp.array([self.xsize,0])
        if L is None:
            self.L = 0.5 * (self.one.type.L + self.two.type.L)
        else:
            self.L = L
        if k is None:
            self.k = 0.5 * (self.one.type.k + self.two.type.k)
        else:
            self.k = k

        self._cached_disp = sp.array([0.,0.])
        self._cached_force = sp.array([0.,0.])

    def __repr__(self):
        return "<Link: k={0.k}, L={0.L}, betwixt:{0.one},{0.two}>".format(self)

    @property
    def calculation_necessary(self):
        if sp.allclose( self.disp, self._cached_disp ):
            return False
        return True

    @property
    def disp(self):
        """ Displacement from one to two """
        direct = self.one.pos - self.two.pos
        around = self.one.pos + self.offset - self.two.pos
        if norm(direct) <= norm(around):
            return direct
        else:
            return around

    @property
    def extension(self):
        """ Get the current extension of the link """
        return norm(self.disp)

    @property
    def energy(self):
        """ Get the energy stored in a link """
        ext = self.extension
        return 0.5*self.k*pow(ext-self.L,2)

    @property
    def force(self):
        """ Get the force the link enacts """
        if self.calculation_necessary:
            ext = self.extension
            disp = self.disp
            self._cached_disp = disp
            force = -self.k * ( ext - self.L ) * unitize(disp)
            self._cached_force = force
            return force
        else:
            return self._cached_force






class Links(object):
    """ Container for Links
        
        Main Attributes:
             * data : holds links, a dictionary where for a pair of cells holds a pointer to the Link for those two cells
             * neighbors : holds neighbor information, a dictionary that for each cell stores a set of its neighbors.
             
        This class basically a wrapper for the builtin dictionary, such that
        accessing its arguments is independent of order.    
    """
    def __init__(self):
        self.data = {}  #where links go. 
        self.neighbors = defaultdict(set) 

    def __repr__(self):
        return "<Links: has {} links betwixt {} cells>".format(len(self.data),
                                                            len(self.neighbors))
    def ord(self,one,two):
        return tuple(sorted((one,two)))

    def add_link(self,one,two,*args,**kwargs):
        """ Add a link between cells one and two """
        self.data[self.ord(one,two)] = Link(one,two,*args,**kwargs)
        self.neighbors[one].add(two)
        self.neighbors[two].add(one)

    def remove_link(self,one,two):
        """ Remove a link between cells one and two """
        del self.data[self.ord(one,two)]
        self.neighbors[one].remove(two)
        self.neighbors[two].remove(one)

    def remove_cell(self,cell):
        """ Remove all references to a cell, both links and neighbors """
        for neigh in self.data[cell]:
            del self.data[self.ord(cell,neigh)]
            self.neighbors[neigh].remove(cell)
        del self.neighbors[cell]

    def get_link(self,one,two):
        """ Get the link between cells one and two, order independent """
        return self.data[self.ord(one,two)]
        
    def get_neighbors(self,cell):
        """ Get the neighbors of cell cell """
        return self.neighbors[cell]

    def iteritems(self):
        return self.data.iteritems()
    def __iter__(self):
        return iter(self.data.values())

    def __getitem__(self,elem):
        if hasattr(elem,'__iter__'):
            one,two = elem
            return self.get_link(one,two)
        else:
            if elem in self.neighbors:
                return self.get_neighbors(elem)
            else:
                raise KeyError

    def __delitem__(self,elem):
        if hasattr(elem,'__iter__'):
            one,two = elem
            self.remove_link(one,two)
        else:
            if elem in self.neighbors:
                self.remove_cell(elem)
            else:
                raise KeyError
    #try to do __getitem__, __setitem__, __delitem__


########################################################
### Simulation Class ###################################
########################################################


class CancerSim:
    """ 
        The main Cancer Simulation Class.
        
        Creates an array of Cells, allows for the designation of cancer cells
        And the evolution of the cells thereafter.
    """
    def __init__(self,boxsize=(XSIZE,YSIZE),basal_height = 10.,seed=None,xi=2.,a=0.1):
        """ Initialize the simulation """

        #set seed
        if seed is None:
            self.seed = int(time.time())
        else:
            self.seed = seed

        sp.random.seed(self.seed)
        random.seed(self.seed)
        #input params
        self.boxsize = boxsize
        self.xi = xi
        self.a = a
        self.XSIZE,self.YSIZE = boxsize
        self.basal_height = basal_height

        #KDTree
        self._kdtree = None

        self._updated = True
        self.T = 0

        # cell types (should be arguments)
        self.cancer = CellType('Cancer',0,0.01,1.0,color='y')
        self.epidermal = CellType('Epidermal',1,0.5,1.0,color='b')
        self.basal = CellType('Basal',2,1.0,0.2,color='k')
        self.dermal = CellType('Dermal',3,0.1,2.0,color='r')

        self.num_cells = 0

        # containers
        self.links = Links()
        self._cell_arr = sp.array([])

        self.cells = []
        self._ghosts = []
        self._ghost_cutoff = 4
        self._ghost_offset = sp.array([boxsize[0],0.])
        self.cancer_cells = []




    def _triang_lattice(self):
        """ Create a triangular grid of points """
        XSIZE, YSIZE = self.boxsize

        #setup the epicells
        epispacing = self.epidermal.L
        xspace,yspace = epispacing , epispacing * sp.sqrt(3)
        for i in sp.arange(0,XSIZE,xspace):
            for ind,j in enumerate(sp.arange(self.basal_height,YSIZE,yspace)):
                if ind:
                    cell1 = Cell([i,j],self.epidermal,self.num_cells)
                    self.add_cell(cell1)
                cell2 = Cell([i+0.5*xspace,j+0.5*yspace],self.epidermal,self.num_cells)
                self.add_cell(cell2)

                #add ghosts for first few layers
                if i<self._ghost_cutoff:
                    if ind:
                        ghost1 = GhostCell(cell1,1,XSIZE)
                        self._ghosts.append(ghost1)
                    ghost2 = GhostCell(cell2,1,XSIZE)
                    self._ghosts.append(ghost2)

                #add ghosts for last few layers
                if i>(XSIZE-self._ghost_cutoff):
                    if ind:
                        ghost1 = GhostCell(cell1,-1,XSIZE)
                        self._ghosts.append(ghost1)
                    ghost2 = GhostCell(cell2,-1,XSIZE)
                    self._ghosts.append(ghost2)

        #setup the bottom cells
        dermalspacing = self.dermal.L
        xspace,yspace = dermalspacing , dermalspacing*sp.sqrt(3)
        for i in sp.arange(0,XSIZE,xspace):
            for j in sp.arange(0,self.basal_height*0.9,yspace):
                cell1 = Cell([i,j],self.dermal,self.num_cells)
                self.add_cell(cell1)
                if j+0.5*yspace < self.basal_height:
                    cell2 = Cell([i+0.5*xspace,j+0.5*yspace],self.dermal,self.num_cells)
                    self.add_cell(cell2)

                #add ghosts for first few layers
                if i<self._ghost_cutoff:
                    ghost1 = GhostCell(cell1,1.,XSIZE)
                    ghost2 = GhostCell(cell2,1,XSIZE)
                    self._ghosts.extend([ghost1,ghost2])

                #add ghosts for last few layers
                if i>(XSIZE-self._ghost_cutoff):
                    ghost1 = GhostCell(cell1,-1,XSIZE)
                    ghost2 = GhostCell(cell2,-1,XSIZE)
                    self._ghosts.extend([ghost1,ghost2])

        #setup the middle cells
        basalspacing = self.basal.L
        for i in sp.arange(0,XSIZE,basalspacing):
            cell = Cell([i,self.basal_height],self.basal,self.num_cells)
            self.add_cell(cell)
            if i<self._ghost_cutoff:
                ghost = GhostCell(cell,1,XSIZE)
                self._ghosts.append(ghost)
            if i>(XSIZE-self._ghost_cutoff):
                ghost = GhostCell(cell,-1,XSIZE)
                self._ghosts.append(ghost)


    
    def get_pos_arr(self,force=False):
        """ Get an array of all of the cell positions """
        #if self._updated is False or force:
        #    return self._cell_arr
        
        self._cell_arr = sp.zeros((len(self.cells),2))
        for (i,cell) in enumerate(self.cells):
            self._cell_arr[i] = cell.pos
        
        self._updated = False
        return self._cell_arr

    def get_radius_arr(self):
        rad_arr=sp.zeros(len(self.cells))
        for (i,cell) in enumerate(self.cells):
            rad_arr[i] = cell.radius
        return rad_arr

    def _get_kdtree(self,force=False):
        """ Generate a KDTree for the cells, 
            allows for efficient geometric neighbor computation """
        pos = self.get_pos_arr(force).copy()
        self._kdtree = KDTree(pos)

        return self._kdtree

    def _query_point(self,x,r,eps=None):
        """ Get all of the cell inds near point, with radius r """        
        kdtree = self._get_kdtree()

        if eps:
            cell_inds = kdtree.query_ball_point(x,r,eps)
        else:
            cell_inds = kdtree.query_ball_point(x,r)
        cells = [ self.cells[ind] for ind in cell_inds ]
        return cells

    def _get_vel_arr(self):
        """ Get an array of all of the cell velocities """
        vel_arr = sp.zeros((self.num_cells,2))
        for (i,cell) in enumerate(self.cells):
            vel_arr[i] = cell.vel
        return vel_arr

    def _update_pos(self,pos_arr):
        """ Update all of the cell positions with an array """
        for (pos,cell) in zip(pos_arr,self.cells):
            cell.pos = pos
        self._cell_arr = pos_arr

    def _update_vel(self,vel_arr):
        """ Update all of the cell velocities with an array """
        for (vel,cell) in zip(vel_arr,self.cells):
            cell.vel = vel
        
    def _get_ghost_pos_arr(self):
        """ Get all of the ghost positions """
        arr = sp.zeros((len(self._ghosts),2))
        for ind,cell in enumerate(self._ghosts):
            arr[ind] = cell.pos
        return arr

    def _update_ghosts(self):
        """ Update the positions of all of the ghost cells """
        for ghost in self._ghosts:
            ghost.update()

    def jiggle(self,sigma=0.1,ghosts=True):
        """ Jiggle the atom positions """
        pos = self.get_pos_arr()

        sigarr = sp.array([cell.type.L for cell in self.cells])
        randn = sp.randn(self.num_cells,2)

        newpos = pos + sigma*(sigarr*randn.T).T

        self._update_pos(newpos)

        if ghosts:
            self._update_ghosts()

    def _set_radii(self):
        """ set radii as the average of the links starting from each cell """
        for cell in [cell for cell in self.cells if cell.type == self.epidermal]:
            average_length=0.0
            norm=0.0
            for neigh in self.links.get_neighbors(cell):
                average_length += self.links.get_link(cell,neigh).L/2.0
                norm += 1.

            cell.radius=average_length/norm
            
        for cell in [cell for cell in self.cells if cell.type == self.dermal]:
             cell.radius=self.epidermal.L/2.0

    def _set_radii_min(self):
        """ set radii as the average of the links starting from each cell """
        for cell in [cell for cell in self.cells if cell.type == self.epidermal]:
            min_length=100000.
            for neigh in self.links.get_neighbors(cell):
                if min_length > self.links.get_link(cell,neigh).L/2.0:
                    min_length = self.links.get_link(cell,neigh).L/2.0
            
            cell.radius=min_length
            
        for cell in [cell for cell in self.cells if cell.type == self.dermal]:
             cell.radius=self.epidermal.L/2.0

    def _freeze_links(self):
        """ Adjust all of the links to be their current extension """
        for link in self.links:
            link.L = link.extension

        self._set_radii()


    def _filter_ghosts(self,one,two):
        if isinstance(one,GhostCell) and isinstance(two,GhostCell):
            raise Exception("DoubleGhost")
        elif isinstance(one,GhostCell):
            return one.original,two
        elif isinstance(two,GhostCell):
            return one,two.original
        else:
            return one,two


    def _clear_links(self):
        """ Clear all Links """
        self.links = Links()

    def delaunay(self):
        """ Delaunay routine, sets the initial links """

        #first get the positions of all the cells and the ghosts
        num_cells = len(self.cells)
        num_ghosts = len(self._ghosts)
        fulllist = self.cells + self._ghosts
        num_full = len(fulllist)

        arr = sp.zeros((num_full,2))
        for ind,cell in enumerate(fulllist):
            arr[ind] = cell.pos
        
        #get the Delaunay construction
        tri = Delaunay(arr)

        #add the links
        for i,j,k in tri.vertices:
            cellone = fulllist[i]
            celltwo = fulllist[j]
            cellthree = fulllist[k]
            try:
                one,two = self._filter_ghosts(cellone,celltwo)
                self.add_bond(one,two)
            except Exception:
                pass
            try:
                one,two = self._filter_ghosts(celltwo,cellthree)
                self.add_bond(one,two)
            except Exception:
                pass
            try:
                one,two = self._filter_ghosts(cellthree,cellone)
                self.add_bond(one,two)
            except Exception:
                pass
        
    def add_cell(self,cell):
        """ Add the cell: cell """
        self.cells.append(cell)
        self.num_cells += 1

    def add_bond(self,one,two):
        """ Add a bond between cells one and two """
        self.links.add_link(one,two,xsize=self.XSIZE)

    def remove_bond(self,one,two):
        """ Remove a bond between cells one and two """
        self.links.remove_link(one,two)

    def remove_cell(self,cell):
        """ Remove the cell: cell, and all bonds for that cell """
        self.cells.remove(cell)
        self.links.remove_cell(cell)

    def get_neighbors(self,cell):
        """ Get the linked neighbor cells of cell """
        return self.links.get_neighbors(cell)

    def add_cancer_cell(self,x,r,eps=None):
        """ randomly make a cell a cancer cell """
        cells = self._query_point(x,r,eps)
        if cells:
            cell = random.choice(cells)
            self.cancer_cells.append(cell)
            cell.type = self.cancer
            try:
                self.links.remove_cell(cell)
            except KeyError:
                pass
        else:
            raise Exception("No targets found at {} within radius {}".format(x,r))


    def duplicate_cancer_cell(self,cancer=None,disp_frac = 0.15):
        """ Duplicate the cancer cell: cancer """
        if cancer is None:
            cancer = random.choice(self.cancer_cells)
        
        #need to choose a random direction and do the relaxation
        L = disp_frac * cancer.type.L
        theta = sp.rand()*2*sp.pi

        disp = L * sp.array([sp.sin(theta),sp.cos(theta)])

        newcell = Cell(cancer.pos + disp,self.cancer,self.num_cells)
        newcell.radius = cancer.radius

        cancer.pos = cancer.pos - disp

        self.cancer_cells.append(newcell)
        self.add_cell(newcell)

        """
        neighs = self.links.get_neighbors(cancer).copy()

        
        for neigh in neighs:
            link_disp =  neigh.pos - cancer.pos
            if sp.vdot(link_disp,disp) >= 0:
                #remove old link, create new one.
                self.links.remove_link(cancer,neigh)
                self.links.add_link(newcell,neigh)
        """

        #self.links.add_link(newcell,cancer)

        self._updated = True


    def time_step(self):
        """ Run a time step, duplicate a cancer cell,
             do a FIRE relaxation, and plot """
        self.duplicate_cancer_cell()
        self.fire()
        self.plot_cells()


    def plot_cells(self,clf=True,fignum=1,ghosts=False,*args,**kwargs):
        """ Plot the current configuration """
        pos_arr = self.get_pos_arr()

        py.figure(fignum)
        if clf:
            py.clf()
             
        py.scatter(pos_arr[:,0],pos_arr[:,1],
                        c=[i.type.color for i in self.cells],
                        s=50,
                        zorder=10,
                        *args,**kwargs)
        if ghosts:
            ghost_arr = self._get_ghost_pos_arr()
            py.scatter(ghost_arr[:,0],ghost_arr[:,1],
                        c = [i.original.type.color for i in self._ghosts],
                        s = 30,
                        zorder=10,
                        alpha = 0.3,
                        *args,**kwargs)
        py.axis('equal')

    def my_circle_scatter(self, axes, x_array, y_array, rad_array, col_array, **kwargs):
           for x, y, R, c in zip(x_array, y_array , rad_array, col_array):
               circle = py.Circle((x,y), radius=R, color = c, **kwargs)
               axes.add_patch(circle)
           return True



    def plot_sized_cells(self,clf=True,fignum=1,ghosts=False,*args, **kwargs):
        """ Plot the current configuration using circles"""
        pos_arr = self.get_pos_arr()
        rad_arr = self.get_radius_arr()
        col_arr = [i.type.color for i in self.cells]

        py.figure(fignum)
        if clf:
            py.clf()
           
        axes=py.axes()
        self.my_circle_scatter(axes,pos_arr[:,0],pos_arr[:,1],rad_arr, col_arr,**kwargs)
                     
        if ghosts:
            ghost_arr = self._get_ghost_pos_arr()
            py.scatter(ghost_arr[:,0],ghost_arr[:,1],
                        c = [i.original.type.color for i in self._ghosts],
                        s = 30,
                        zorder=10,
                        alpha = 0.3,
                        *args,**kwargs)
        py.axis('equal')


    def plot_links(self,clf=False,cutoff=None,fignum=1,ghosts=False,*args,**kwargs):
        """ Plot the links between cells """  
        if cutoff is None:
            cutoff = self.XSIZE/2.          
        py.figure(fignum)
        if clf:
            py.clf()
        
        for link in self.links:
            d12=link.one.pos-link.two.pos
            abs_d12=norm(d12)
            if abs_d12 < cutoff:
                data = sp.array([ link.one.pos, link.two.pos ])
                py.plot(data[:,0],data[:,1],
                            c=py.cm.jet( min(link.energy*30.,1.) ),
                            alpha=0.6,
                            *args, **kwargs )

    @property
    def forces(self):
        """ get the forces between cells, as array, both from links
            and from the native force_func
        """
        pos = self.get_pos_arr(force=True)

        force_arr = sp.zeros_like(pos)

        for link in self.links:
            force = link.force
            force_arr[link.one.index] += force
            force_arr[link.two.index] -= force


        kdtree = self._get_kdtree(force=True)
        for i,j in kdtree.query_pairs(self.xi*1.0):
            
            force = self.force_func_hertz(self.cells[i], self.cells[j] )
            #disp = self.cells[i].pos - self.cells[j].pos
            #L = norm(disp)
            #force = 2 * self.a**4 * ( 2 * self.xi**2 - 3 * self.xi * L + L**2 )/( self.xi**2 * L**6 ) * disp
            force_arr[i] += force
            force_arr[j] -= force

        return sp.nan_to_num(force_arr)


    def force_func(self,cell1,cell2):
        """ the native force function between two positions """
        x1 = cell1.pos
        x2 = cell2.pos
        disp = x1 - x2
        mod_disp = norm(disp)
        force = 2 * self.a**4 * ( 2 * self.xi**2 - 3 * self.xi * mod_disp + mod_disp**2 )/( self.xi**2 * mod_disp**6 ) * disp
        
        return force

    def force_func2(self,cell1,cell2):
        """ the native force function between two positions """
        x1 = cell1.pos
        x2 = cell2.pos
        r1 = cell1.radius
        r2 = cell2.radius
        disp = x1 - x2
        mod_disp = norm(disp)
        a1=self.a*(r1+r2)
        xi1=self.xi*(r1+r2)
        force = 2 * a1**4 * ( 2 * xi1**2 - 3 * xi1 * mod_disp + mod_disp**2 )/( xi1**2 * mod_disp**6 ) * disp
        
        return force

    def force_func_hertz(self,cell1,cell2):
         """ the Hertz force between two cells """
         x1 = cell1.pos
         x2 = cell2.pos
         r1 = cell1.radius
         r2 = cell2.radius
         disp = x1 - x2
         mod_disp = norm(disp)
         delta=(r1+r2)-mod_disp
         if delta > 0.0:
             force = self.a*delta**1.5*disp/mod_disp
         else:
             force= 0.0

         return force


    @property
    def energy(self):
        """ get the energy of the current configuration """
        tot_energy = 0
        for link in self.links:
            tot_energy += link.energy
        return tot_energy

    def fire(self,fmax=0.1,
            Nmin=5.,finc=1.1,fdec=0.5,alphastart=0.1,fa=0.99,deltatmax=10.,
            maxsteps = 10**5):
        """ Do a fire relaxation """

        alpha = alphastart
        deltat = 0.1

        pos = self.get_pos_arr(force=True)
        v = sp.zeros_like(pos)
        self._update_vel(v)

        v = self._get_vel_arr()

        steps_since_negative = 0

        def norm_arr(vec):
            return sp.sqrt(sp.sum(vec**2,1))
        def unitize_arr(vec):
            return ((vec.T)/norm(vec)).T

        forces = sp.nan_to_num(sp.array([ [sp.inf,sp.inf]]))

        step_num = 0

        print "Beginning FIRE Relaxation -- fmax={}".format(fmax)

        while max(norm_arr(forces)) > fmax and step_num < maxsteps:
            forces = self.forces

            power = sp.vdot(forces,v)
            print "Step: {}, max_force: {}, power: {}".format(step_num,max(norm_arr(forces)), power)

            v = (1.0 - alpha)*v + alpha*(norm_arr(v)*unitize_arr(forces).T).T

            if power>0.:
                if steps_since_negative > Nmin:
                    deltat = min(deltat * finc, deltatmax)
                    alpha = alpha*fa
                steps_since_negative += 1

            else:
                steps_since_negative = 0

                deltat = deltat * fdec
                v *= 0.
                alpha = alphastart

            v += forces*deltat
            pos += v*deltat
            self._update_pos(pos)
            step_num += 1
        
        self._update_pos(pos)
        self._update_vel(v)
        print "Relaxation finished..."
    
if __name__ == "__main__":
    Q = CancerSim()
    Q._triang_lattice()
    Q.delaunay()
    Q._freeze_links()

    Q.add_cancer_cell([XSIZE/2.,YSIZE/2 + 3],1)

    Q.plot_cells()


    self = Q


"""
TODO:  have links know about periodic boundary conditions (maybe)
        freeze links (DONE)
        Ghost cells need update method.  (DONE)
        fire relaxation (DONE)
        set and divide cancer cells (DONE)
        long range forces (DONE)
            
        cache the link calcs
        cache the KDTree calcs?
        allow more transparent custimization
        expose CellTypes
        use logging module
"""
