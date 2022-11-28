import numpy as np
from ._CFunctions import _CGetCon2020Params,_CSetCon2020Params
import ctypes as ct

def _GetCFG():
	'''
	Get the current config dictionary
	
	
	'''
	eqtype = ct.c_char_p("        ".encode('utf-8'))
	mui = np.zeros(1,dtype='float64')
	irho = np.zeros(1,dtype='float64')
	r0 = np.zeros(1,dtype='float64')
	r1 = np.zeros(1,dtype='float64')
	d = np.zeros(1,dtype='float64')
	xt = np.zeros(1,dtype='float64')
	xp = np.zeros(1,dtype='float64')
	Edwards = np.zeros(1,dtype='bool')
	ErrChk = np.zeros(1,dtype='bool')
	CartIn = np.zeros(1,dtype='bool')
	CartOut = np.zeros(1,dtype='bool')
	Smooth = np.zeros(1,dtype='bool')
	DeltaRho = np.zeros(1,dtype='float64')
	DeltaZ = np.zeros(1,dtype='float64')

	g = np.zeros(1,dtype='float64')
	azfunc = ct.c_char_p("          ".encode('utf-8'))
	wO_open = np.zeros(1,dtype='float64')
	wO_om = np.zeros(1,dtype='float64')
	thetamm = np.zeros(1,dtype='float64')
	dthetamm = np.zeros(1,dtype='float64')
	thetaoc = np.zeros(1,dtype='float64')
	dthetaoc = np.zeros(1,dtype='float64')


	_CGetCon2020Params(mui,irho,r0,r1,d,xt,xp,eqtype,Edwards,ErrChk,
						CartIn,CartOut,Smooth,DeltaRho,DeltaZ,g,azfunc,
						wO_open,wO_om,thetamm,dthetamm,thetaoc,dthetaoc)
	
	cfg = {}
	cfg['mu_i'] = mui[0]
	cfg['i_rho'] = irho[0]
	cfg['r0'] = r0[0]
	cfg['r1'] = r1[0]
	cfg['d'] = d[0]
	cfg['xt'] = xt[0]
	cfg['xp'] = xp[0]
	cfg['Edwards'] = Edwards[0]
	cfg['error_check'] = ErrChk[0]
	cfg['CartesianIn'] = CartIn[0]
	cfg['CartesianOut'] = CartOut[0]
	cfg['equation_type'] = eqtype.value.decode()
	cfg['Smooth'] = Smooth[0]
	cfg['DeltaRho'] = DeltaRho[0]
	cfg['DeltaZ'] = DeltaZ[0]

	cfg['g'] = g[0]
	cfg['azfunc'] = azfunc.value.decode()
	cfg['wO_open'] = wO_open[0]
	cfg['wO_om'] = wO_om[0]
	cfg['thetamm'] = thetamm[0]
	cfg['dthetamm'] = dthetamm[0]
	cfg['thetaoc'] = thetaoc[0]
	cfg['dthetaoc'] = dthetaoc[0]

	return cfg

def _SetCFG(cfg):
	'''
	Set the model config using a dictionary.
	
	'''
	
	eqtype = ct.c_char_p(cfg['equation_type'].encode('utf-8'))
	mui = np.array([cfg['mu_i']],dtype='float64')
	irho = np.array([cfg['i_rho']],dtype='float64')
	r0 = np.array([cfg['r0']],dtype='float64')
	r1 = np.array([cfg['r1']],dtype='float64')
	d = np.array([cfg['d']],dtype='float64')
	xt = np.array([cfg['xt']],dtype='float64')
	xp = np.array([cfg['xp']],dtype='float64')
	Edwards = np.array([cfg['Edwards']],dtype='bool')
	ErrChk = np.array([cfg['error_check']],dtype='bool')
	CartIn = np.array([cfg['CartesianIn']],dtype='bool')
	CartOut = np.array([cfg['CartesianOut']],dtype='bool')
	Smooth = np.array([cfg['Smooth']],dtype='bool')
	DeltaRho = np.array([cfg['DeltaRho']],dtype='float64')
	DeltaZ = np.array([cfg['DeltaZ']],dtype='float64')

	g = np.array([cfg['g']],dtype='float64')
	azfunc = ct.c_char_p(cfg['azfunc'].encode('utf-8'))
	wO_open = np.array([cfg['wO_open']],dtype='float64')
	wO_om = np.array([cfg['wO_om']],dtype='float64')
	thetamm = np.array([cfg['thetamm']],dtype='float64')
	dthetamm = np.array([cfg['dthetamm']],dtype='float64')
	thetaoc = np.array([cfg['thetaoc']],dtype='float64')
	dthetaoc = np.array([cfg['dthetaoc']],dtype='float64')



	
	_CSetCon2020Params(mui,irho,r0,r1,d,xt,xp,eqtype,Edwards,ErrChk,
						CartIn,CartOut,Smooth,DeltaRho,DeltaZ,g,azfunc,
						wO_open,wO_om,thetamm,dthetamm,thetaoc,dthetaoc)
						
						
def Config(*args,**kwargs):
	'''
	Set and return the Con2020 model configuration
	
	Input Arguments
	===============
	The only accepted argument here is the string "default", i.e.
	Config("default") which sets the model to default settings.
	
	Keywords (shorthand keywords in brackets)
	=========================================
	mu_i_div2__current_parameter_nT (mu_i): float
		mu0i0/2 term (azimuthal current sheet parameter), in nT
	i_rho__radial_current_MA (i_rho) : float
		radial current term from Connerney et al., 2020
		NOTE: The default value (16.7 MA) is the average value from
		Connerney et al 2020. This value was shown to vary from one 
		pass to the next, where Table 2 provides radial current 
		density values for 23 of the first 24
		perijoves. (Only used for Connerney model for the azimuthal field)
	r0__inner_rj (r0) : float
		Inner edge of current disk in Rj
	r1__outer_rj (r1) : float
		Outer edge of current disk in Rj
	d__cs_half_thickness_rj (d) : float
		Current sheet half thickness in Rj
	xt__cs_tilt_degs (xt) : float
		Current sheet tilt in degrees
	xp__cs_rhs_azimuthal_angle_of_tilt_degs (xp) : float
		Current sheet tilt longitude (right handed) in degrees
	Edwards : bool
		Use Edwards et al equations.
	equation_type: str
		Define method for calculating the current sheet field, may be 
		one of the following: 'hybrid'|'analytic'|'integral'
		See notes below for more information.
	error_check : bool
		If True (default) then inputs will be checked for potential errors.		
	CartesianIn : bool
		If True (default) the inputs to the model will be expected to be 
		in Cartesian right-handed System III coordinates. If False, then
		the inputs should be in spherical polar coordinates.
	CartesianOut : bool
		If True (default) the output magnetic field will be in Cartesian
		right-handed System III coordinates. Otherwise, the magnetic 
		field components produced will be radial, meridional and 
		azimuthal.	
	Smooth : bool
		If True, then smooth transitions are used about r0 and r1
		(currently Edwards only).
	g : float
		Dipole coefficient, nT. (Used for L-MIC model)
	azfunc : str
		"connerney" - uses the Connerney et al 2020 method for calculating
					azimuthal component of the magnetic field.
		or
		"lmic" - uses the Leicester Magnetosphere-Ionosphere Coupling model
					(e.g. Cowley et al 2005, 2008) to provide the azimuthal
					component of the magnetic field.
	wO_open	: float
			angular velocity ratio of open flux to planetary spin
	wO_om : float
		angular velocity ratio of outer magnetosphere to planetary spin
	thetamm : float
		ionospheric footprint latitude of the middle magnetosphere (where plasma 
		goes from rigid corotation to subcorotation) in radians.
	dthetamm : float
		width of the middle magnetosphere in radians.
	thetaoc	: float
		ionospheric latitude of the open-closed field line boundary, in radians.
	dthetaoc : float
		width of the open-closed field line boundary, in radians.
	
	'''

	#list the default arguments here
	defargs = {	'mu_i'			: 139.6,
				'i_rho' 		: 16.7,
				'r0'			: 7.8,
				'r1'			: 51.4,
				'd'				: 3.6,
				'xt'			: 9.3,
				'xp'			: 155.8,
				'equation_type'	: 'hybrid',
				'Edwards'		: True,
				'error_check'	: True,
				'CartesianIn'	: True,
				'CartesianOut'	: True,
				'Smooth'		: False,
				'DeltaRho'		: 1.0,
				'DeltaZ'		: 0.1,
				'g'				: 417659.3836476442,
				'azfunc'		: "connerney",
				'wO_open'		: 0.1,
				'wO_oc'			: 0.35,
				'thetamm'		: 16.1,
				'dthetamm'		: 0.5,
				'thetaoc'		: 10.716,
				'dthetaoc'		: 0.125,}
				
	
	if len(args) == 1:
		if args[0] == "default":
			#set the configuration of the model to default
			_SetCFG(defargs)
			#note that we can still modify the default configuration
			#using the keywords provided after doing this		
	
	#return the current configuration
	cfg = _GetCFG()

	#list the long names
	longnames = {	'mu_i'	: 'mu_i_div2__current_parameter_nT',
					'r0'	: 'r0__inner_rj',
					'r1'	: 'r1__outer_rj',
					'd'		: 'd__cs_half_thickness_rj',
					'xt'	: 'xt__cs_tilt_degs',
					'xp'	: 'xp__cs_rhs_azimuthal_angle_of_tilt_degs',
					'i_rho'	: 'i_rho__radial_current_MA'		  }
						
	#check input kwargs
	#for those which exist (either in long or short name form) add
	#them to this object using the short name as the object tag
	#Otherwise use the existing value
		
	#the input keys
	ikeys = list(kwargs.keys())
		
	#current config keys
	ckeys = list(cfg.keys())
		
	#short and long name keys
	skeys = list(longnames.keys())
	lkeys = [longnames[k] for k in skeys]

		
	#loop through each one		
	for k in ckeys:
		if k in ikeys:
			#short name found in kwargs - add to this object
			cfg[k] = kwargs[k]
		elif longnames.get(k,'') in ikeys:
			#long name found - add to object
			cfg[k] = kwargs[longnames[k]]
		else:
			#key not found, use existing parameter
			pass
		
	#check for additional keys and issue a warning
	for k in ikeys:
		if not ((k in skeys) or (k in lkeys) or (k in ckeys)):
			print("Keyword argument {:s} unrecognized, ignoring.".format(k))

	#update the configuration
	_SetCFG(cfg)
	
	#get a copy of the config stored in C++
	cfg = _GetCFG()

	
	return cfg
