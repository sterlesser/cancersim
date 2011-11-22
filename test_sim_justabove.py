from code import cancer
import pylab as py

from config import config

XSIZE = config['XSIZE']
YSIZE = config['YSIZE']
config['first_cancer_cell_yoffset']=1


py.figure(1)
#py.show()

Q = cancer.CancerSim(config)
Q._setup()

py.ioff()

Q.plot_sized_cells()
py.draw()

imagepath = "justabove/{:06d}.png"


for i in range(1000):
    Q.time_step()
    #Q.plot_links()
    py.savefig(imagepath.format(i))
    Q.save('justabove/state.dat')
    #Q.plot_links()
    #py.xlim((5,16))
    #py.ylim((9,17))
    #py.draw()
      
    

