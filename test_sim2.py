from code import cancer
import pylab as py

XSIZE = 10
YSIZE = 10

py.figure(1)
py.show()

Q = cancer.CancerSim(boxsize=[XSIZE,YSIZE],basal_height = 4.5,a=0.5,xi=0.5,seed=None)
Q._triang_lattice()
Q.jiggle(sigma=0.1)
Q.delaunay()
Q._freeze_links()

Q.add_cancer_cell([XSIZE/2.,YSIZE/2 + 0.5],1)

Q.plot_cells()
py.draw()


for i in range(50):
    Q.time_step()
    #Q.plot_links()
    #py.xlim((5,16))
    #py.ylim((9,17))
    py.draw()
      
    

