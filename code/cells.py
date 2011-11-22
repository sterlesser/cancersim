
import scipy as sp 

from logger import logger
import pprint

base_logger = logger.getChild('cells')
base_logger.info('Inside cells.py')
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
    
    logger = base_logger.getChild('CellType')
    def __init__(self,name,type_ind,k,L,color='b',maxstretch=1.2,*args,**kwargs):
        self.name = name
        self.type_ind = type_ind 
        self.k = k 
        self.L = L  
        self.color = color
        self.maxstretch = maxstretch
        
        logger.debug("""Creating a CellType instance with:
                    {info}""".format(info=pprint.pformat(self.__dict__)))

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
    logger = base_logger.getChild('Cell')

    def __init__(self,pos,cell_type,index=0):
        self.pos = sp.array(pos)
        self.type = cell_type
        self.vel = sp.zeros_like(self.pos)
        self.index = index
        self.radius = self.type.L/2.0
        
        logger.debug("""Creating a Cell with:
                            {info}""".format(info=pprint.pformat(self.__dict__)))

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

    logger = base_logger.getChild('GhostCell')
    def __init__(self,original,xsize,sign=1):
        self.offset = sign* sp.array([xsize,0.])
        self.sign = sign
        self.original = original
        self.type = original.type
        self.vel = original.vel

        self.pos = self.original.pos + self.offset

        logger.debug("""Creating a GhostCell with:
                        {info}""".format(info=pprint.pformat(self.__dict__)))
                            

    def update(self):
        """ Update the position of the GhostCell based on the offset and sign """
        self.pos = self.original.pos + self.offset

    def __repr__(self):
        return "<GhostCell: pos:{0.pos} of cell {0.original}>".format(self)
