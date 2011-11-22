from code import cancer
import pylab as py

from config import config

config['XSIZE'] = 10
config['YSIZE'] = 10
config['basal_height'] = 6
config['first_cancer_cell_yoffset'] = 1.5


py.figure(1)
py.show()

Q = cancer.CancerSim(config)
Q._setup()

#py.ioff()

Q.plot_sized_cells()
py.draw()



def run(steps=10):
	for i in range(steps):
	    Q.time_step()
	    #Q.plot_links()
	    py.xlim((0,config['XSIZE']))
	    py.ylim((0,config['YSIZE']))
	    py.draw()
      
if __name__ == "__main__":
	run(10)

