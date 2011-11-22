ToDo
=====

Force function that depends on the type of the cell

Cell packing, long links at top and bottom.
    - could push from the top?
    - look for boundary condition settings in Delaunay - Qhull
    - could cut outliers

Take a look at the boundary conditions, make sure they're working.

Volume Preserving
    -http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1001045
    
Pass in Cell Types

Dictionary configuration?

Run a large sample

Exchange springs for adhesion
    - turns on when they get close, has finite range
    
Epidermis is packed structure, dermis is spring links,
     basal cells, just have short range repulsions, not cells but nodes of 
        a network, if stretched too much breaks


WEEK2
-----------------

Want the Basal membrane to be an impenetrable barrier.  Strong repulsion that is short ranged.

Basal cells don't interact.

Play around with the force constants

Different constants for membrane interactions and other cell interactions

Configuration dictionary.

Convergence issues, maybe don't create cells in center

Break links if they stretch 1.2 times the original length

Figure out what's wrong with the plot sized cell routine

WEEK2 TODO
---------------

Extract parameters.  Look at Plot Sized Cells.  Convergence stuff.  Cell breaking.  

WEEK 3
--------

Look at bio paper to get some physical parameters

Explore FIRE parameters

Try a large run.

Pull link breaking criteria out to config.

Maybe one (or two) link per cell to the Basal membrane

Need to look up some parameters and run some simulations.
Another paper:  
	http://www.cell.com/biophysj/abstract/S0006-3495%2805%2973087-3
	http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1001045
    
WEEK 4
---------

Table 1 of paper 1 has the good parameters.
Check the size of the epidermis.
Cancer cells near the edges split, growth rate according to density
Seed on the base of the membrane
X - plot sized cells??
Voronoi cells paper.
Run the simulation starting in the epidermis and the dermis

