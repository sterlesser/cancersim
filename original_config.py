
config = {
	#General params
	'XSIZE': 10,
	'YSIZE': 10,
	'basal_height': 4.5,
	'seed': None,
	'force_magnitude': 0.1,
	'force_cutoff': 2.,
	'jiggle_sigma': 0.1,

	#First Cancer Cell
	'first_cancer_cell_yoffset': 1.,
	'first_cancer_cell_radius': 1.,

	#CELL PARAMETERS
	'cancer_cell_params': { 
		'name':'Cancer',
		'type_ind':0,
		'k': 0.01,
		'L': 1.,
		'color': 'y'
	},

	'epidermal_cell_params': {
		'name': 'Epidermal',
		'type_ind': 1,
		'k': 0.1,
		'L': 1.0,
		'color': 'b'
	},

	'basal_cell_params': {
		'name': 'Basal',
		'type_ind': 2,
		'k': 0.5,
		'L': 0.2,
		'color': 'k'
	},

	'dermal_cell_params': {
		'name': 'Dermal',
		'type_ind': 3,
		'k': 0.01,
		'L': 2.0,
		'color': 'r'
	},


	#Parameters for FIRE
	'fmax': 0.1,
	'Nmin': 5.,
	'finc': 1.1,
	'fdec': 0.5,
	'alphastart': 0.1,
	'fa': 0.99,
	'deltatmax': 10.,
    'maxsteps': 10**5,
}