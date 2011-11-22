from __future__ import division

""" 
Configuration
this sets the parameters for the simulation


As for units, we are working in units where our natural length scale is 1e-5 m or 10 um

and so that 
A = 2/3( E/(1-nu^2)) = 2/3 ( (1 kPa) / ( 1 - (1/3)^2 ) = 3/4 kPa = 1 

And with a natural adhesion between cells of 200 uN/m, we have the factor below
"""


factor = 200.0/(0.75e4)
config = {
	#General params
	'XSIZE': 40,
	'YSIZE': 30,
	'basal_height': 20,
	'seed': None,
	'force_magnitude': 1.0,  #This is 1/D
	#  This =  2/3 ( E / (1 - nu^2 )  )

	'force_magnitude_basal': 0.10,
	'force_cutoff': 2.0,
	'force_cutoff_basal': 0.4,
	'jiggle_sigma': 0.1,

	#First Cancer Cell
	'first_cancer_cell_yoffset': 1.,
	'first_cancer_cell_radius': 1.,

	#CELL PARAMETERS
	'cancer_cell_params': { 
		'name':'Cancer',
		'type_ind':0,
		'k': 0.*factor,
		'L': 1.,
		'color': 'y',
		'maxstretch': 1.2,
	},

	'epidermal_cell_params': {
		'name': 'Epidermal',
		'type_ind': 1,
		'k': 1.0*factor,
		'L': 1.0,
		'color': 'b',
		'maxstretch': 1.2,
	},

	'basal_cell_params': {
		'name': 'Basal',
		'type_ind': 2,
		'k': 3.0*factor,
		'L': 0.2,
		'color': 'k',
		'maxstretch': 1.5,
	},

	'dermal_cell_params': {
		'name': 'Dermal',
		'type_ind': 3,
		'k': 0.1*factor,
		'L': 2.0,
		'color': 'r',
		'maxstretch': 1.2,
	},


	#Parameters for FIRE
	'fmax': factor/10.,
	'Nmin': 5.,
	'finc': 1.1,
	'fdec': 0.5,
	'alphastart': 0.1,
	'fa': 0.99,
	'deltatmax': 10.,
    'maxsteps': 10**3,
}
