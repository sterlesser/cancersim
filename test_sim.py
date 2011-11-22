from code import cancer
import pylab as py

from config import config

XSIZE = config['XSIZE']
YSIZE = config['YSIZE']

py.figure(1)
py.show()

Q = cancer.CancerSim(config)
Q._setup()

#py.ioff()

Q.plot_sized_cells()
py.draw()

imagepath = "newpics/{:06d}.png"

#py.ioff()

for i in range(1000):
    Q.time_step()
    Q.plot_links()
    py.savefig(imagepath.format(i))
    #Q.plot_links()
    #py.xlim((5,16))
    #py.ylim((9,17))
    py.draw()
      
    

