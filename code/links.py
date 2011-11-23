
from helper import norm, unitize 

from collections import defaultdict
from math import pow

import scipy as sp

import pprint
from logger import logger
base_logger = logger.getChild('links')
base_logger.info('Inside links.py')

########################################################
### Link Stuff  ########################################
########################################################

XSIZE=20

MAXSTRETCH = 1.2



###########
# custom allclose
###########


from scipy import allclose

def allclose2(a,b,rtol=1e-05,atol=1e-08):
    ax,ay = a
    bx,by = b

    return ( abs(ax-bx) <= atol + rtol*abs(bx) ) and ( abs(ay-by) <= atol + rtol*abs(by))



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
    logger = base_logger.getChild('Link')
    def __init__(self,one,two,L=None,k=None,xsize=XSIZE,maxstretch=None):
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
        if maxstretch is None:
            self.maxstretch = 0.5 * (self.one.type.maxstretch + self.two.type.maxstretch)

        self.broken = False

        self._cached_disp = sp.array([0.,0.])
        self._cached_force = sp.array([0.,0.])

        logger.debug("""Created a Link with:
                        {info}""".format(info=pprint.pformat(self.__dict__)))

    def __repr__(self):
        return "<Link: k={0.k}, L={0.L}, betwixt:{0.one},{0.two}>".format(self)

    @property
    def calculation_necessary(self):
        if allclose2( self.disp, self._cached_disp ):
            return False
        return True

    @property
    def disp(self):
        """ Displacement from one to two """
        disp = self.one.pos - self.two.pos
        if norm(disp) > self.xsize/2:
            disp = self.one.pos + self.offset - self.two.pos
        return disp

    @property
    def disp_old(self):
        """ Displacement from one to two """
        direct = self.one.pos - self.two.pos
        around = self.one.pos + self.offset -self.two.pos
        if norm(direct) <= norm(around):
            return direct
        else:
            return around

    def extension_without_breaking(self):
        """ Get the extension of the current link without breaking """
        length = norm(self.disp)

        return length

    @property
    def extension(self):
        """ Get the current extension of the link """
        length = self.extension_without_breaking()
        if not self.broken and length > self.maxstretch * self.L and self.k > 0:
            logger.warning('One of our links is breaking!')
            self.k = 0
            self.broken = True

        return length

    @property
    def energy(self):
        """ Get the energy stored in a link """
        ext = self.extension
        return 0.5*self.k*pow(ext-self.L,2)

    @property
    def force(self):
        """ Get the force the link enacts """
        if self.broken:
            return 0
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
    logger = base_logger.getChild('Links')

    def __init__(self):
        self.data = {}  #where links go. 
        self.neighbors = defaultdict(set) 
        
        logger.debug('Links collection created.')

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

        logger.debug('Created a link between {one} and {two}'.format(one=one,two=two))

    def remove_link(self,one,two):
        """ Remove a link between cells one and two """
        del self.data[self.ord(one,two)]
        self.neighbors[one].remove(two)
        self.neighbors[two].remove(one)

        logger.debug('Removed a link between {one} and {two}'.format(one=one,two=two))

    def remove_cell(self,cell):
        """ Remove all references to a cell, both links and neighbors """
        
        links = self.get_links(cell)
        for link in links:
            self.remove_link(link.one,link.two)

        del self.neighbors[cell]

        logger.debug('Removed the cell {cell}'.format(cell=cell))

    def get_link(self,one,two):
        """ Get the link between cells one and two, order independent """
        return self.data[self.ord(one,two)]
        
    def get_neighbors(self,cell):
        """ Get the neighbors of cell cell """
        return self.neighbors[cell]

    def get_links(self,cell):
        """ Get all of the links that involve cell """
        neighbors = self.get_neighbors(cell)
        links = [self.get_link(cell,neigh) for neigh in neighbors]
        return links
        
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

