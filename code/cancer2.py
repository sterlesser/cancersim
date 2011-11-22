#Cancer Sim

import scipy as sp
import pylab as py

from collections import defaultdict


from scipy.spatial.distance import euclidean
from math import pow


XSIZE = 20
YSIZE = 20

def norm(vec):
    """ Normalize a vector """
    return sp.sqrt(sp.sum(vec**2))

def unitize(vec):
    """ Construct a unit vector from a vector """
    return vec/norm(vec)

class CellType(object):
    """ A holder for some cell type information """
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
        
        Attributes:
            pos: a position array
            type: a type declaration
    """
    def __init__(self,pos,cell_type,index=0):
        self.pos = sp.array(pos)
        self.type = cell_type
        self.vel = sp.zeros_like(self.pos)
        self.index = index

    def __repr__(self):
        return "<Cell type:{0.type.name} pos:{0.pos} ID:{1}>".format(self,id(self))


class GhostCell(Cell):
    """ A ghost cell """
    def __init__(self,original,sign=1,xsize=XSIZE):
        self.offset = sign* sp.array([XSIZE,0.])
        self.sign = sign
        self.original = original
        self.type = original.type
        self.vel = original.vel

        self.pos = self.original.pos + self.offset
    def update(self):
        self.pos = self.original.pos + self.offset

    def __repr__(self):
        return "<GhostCell: pos:{0.pos} of cell {0.original}>".format(self)

class Link(object):
    """ A link object keeps a link between two cells.
    
        Link(cell_one,cell_two,length,k)

        Properties:
            disp - displacement of spring
            extension - compute the current link extension
            energy - compute the current link energy
            force - force of spring
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
        """ Get the force from one to two """
        if self.calculation_necessary:
            ext = self.extension
            disp = self.disp
            self._cached_disp = disp
            force = -self.k * ( ext - self.L ) * unitize(disp)
            self._cached_force = force
            return force
        else:
            return self._cached_force



from collections import defaultdict


class Links(object):
    """ Container for Links, data holds links,
             and neighbors holds neighbors
             
        Methods:
            add_link(one,two,k=None,L=None):
                 adds a link between one and two
                 by default, k and L are average of cell.type vals
            remove_link(one,two): removes links between one and two
            remove_cell(cell): remove a cell and all it's links
            get_link(one,two): get link between one and two
            get_neighbors(cell): get a set of all cells connected to
                cell
            iteritems(): iteritems the links
            __iter__(): iters on the links
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
        self.data[self.ord(one,two)] = Link(one,two,*args,**kwargs)
        self.neighbors[one].add(two)
        self.neighbors[two].add(one)

    def remove_link(self,one,two):
        del self.data[self.ord(one,two)]
        self.neighbors[one].remove(two)
        self.neighbors[two].remove(one)

    def remove_cell(self,cell):
        for neigh in self.data[cell]:
            del self.data[self.ord(cell,neigh)]
        del self.neighbors[cell]

    def get_link(self,one,two):
        return self.data[self.ord(one,two)]
        
    def get_neighbors(self,cell):
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


from scipy.spatial import Delaunay
from scipy.spatial import KDTree
import random
import time

class CancerSim:

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
        self.cancer = CellType('Cancer',0,0.5,1.0,color='y')
        self.epidermal = CellType('Epidermal',1,0.5,1.0,color='b')
        self.basal = CellType('Basal',2,1.0,0.5,color='k')
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
        #if self._updated is False or force:
        #    return self._cell_arr
        
        self._cell_arr = sp.zeros((len(self.cells),2))
        for (i,cell) in enumerate(self.cells):
            self._cell_arr[i] = cell.pos
        
        self._updated = False
        return self._cell_arr

    def _get_kdtree(self,force=False):
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
        vel_arr = sp.zeros((self.num_cells,2))
        for (i,cell) in enumerate(self.cells):
            vel_arr[i] = cell.vel
        return vel_arr

    def _update_pos(self,pos_arr):
        for (pos,cell) in zip(pos_arr,self.cells):
            cell.pos = pos
        self._cell_arr = pos_arr

    def _update_vel(self,vel_arr):
        for (vel,cell) in zip(vel_arr,self.cells):
            cell.vel = vel
        
    def _get_ghost_pos_arr(self):
        arr = sp.zeros((len(self._ghosts),2))
        for ind,cell in enumerate(self._ghosts):
            arr[ind] = cell.pos
        return arr

    def _update_ghosts(self):
        for ghost in self._ghosts:
            ghost.update()

    def jiggle(self,sigma=0.1,ghosts=True):
        """ Jiggle the atom positions """
        pos = self.get_pos_arr()

        sigarr = sp.array([cell.type.L for cell in self.cells])
        sigbool = sp.array([cell.type.name!='Basal' for cell in self.cells])*1
        sigarr = sigarr*sigbool
        print sigbool
        randn = sp.randn(self.num_cells,2)
        
        newpos = pos + sigma*(sigarr*randn.T).T

        self._update_pos(newpos)

        if ghosts:
            self._update_ghosts()

    def _freeze_links(self):
        for link in self.links:
            link.L = link.extension

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
        self.links = Links()

    def delaunay(self):
        """ delaunay routine """

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

    def plot_links(self,clf=False,cutoff=None,fignum=1,ghosts=False,*args,**kwargs):
        """ Plot the links between cells """  
        if cutoff is None:
            cutoff = self.XSIZE/2.          

        py.figure(fignum)
        if clf:
            py.clf()
        
        for link in self.links:
            if link.extension < cutoff:
                data = sp.array([ link.one.pos, link.two.pos ])
                py.plot(data[:,0],data[:,1],
                            c=py.cm.jet( min(link.energy*30.,1.) ),
                            alpha=0.6,
                            *args, **kwargs )
    @property
    def forces(self):
        """ get the forces between cells, as array """
        pos = self.get_pos_arr(force=True)

        force_arr = sp.zeros_like(pos)

        for link in self.links:
            force = link.force
            force_arr[link.one.index] += force
            force_arr[link.two.index] -= force


        kdtree = self._get_kdtree(force=True)
        for i,j in kdtree.query_pairs(self.xi*1.0):
            disp = self.cells[i].pos - self.cells[j].pos
            L = norm(disp)
            force = 2 * self.a**4 * ( 2 * self.xi**2 - 3 * self.xi * L + L**2 )/( self.xi**2 * L**6 ) * disp
            force_arr[i] += force
            force_arr[j] -= force

        return sp.nan_to_num(force_arr)


    def force_func(self,x1,x2):
        """ the native force function between two points """

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

        long range forces (KDTree, but with periodic boundary conditions?)
            just test both a position and its periodic image if it's near the edge (very good)
        cache the link calcs

"""
