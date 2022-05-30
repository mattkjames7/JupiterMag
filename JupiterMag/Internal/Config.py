import numpy as np
from ._CFunctions import _CGetInternalCFG,_CSetInternalCFG
import ctypes as ct

def _GetCFG():
	'''
	Get the current config dictionary
	
	
	'''
	Model = ct.c_char_p("                                ".encode('utf-8'))
	CartIn = np.zeros(1,dtype='bool')
	CartOut = np.zeros(1,dtype='bool')
	Degree = np.zeros(1,dtype='int32')

	_CGetInternalCFG(Model,CartIn,CartOut,Degree)
	
	cfg = {}
	cfg['Model'] = Model.value.decode()
	cfg['CartesianIn'] = CartIn[0]
	cfg['CartesianOut'] = CartOut[0]
	cfg['Degree'] = Degree[0]
	return cfg

def _SetCFG(cfg):
	
	Model = ct.c_char_p(cfg['Model'].lower().encode('utf-8'))
	CartIn = np.array(cfg['CartesianIn'],dtype='bool')
	CartOut = np.array(cfg['CartesianOut'],dtype='bool')
	Degree = np.array(cfg['Degree'],dtype='int32')
	
	_CSetInternalCFG(Model,CartIn,CartOut,Degree)
						
						
def Config(*args,**kwargs):
	'''
	Set and return the internal model configuration
	
	Input Arguments
	===============
	The only accepted argument here is the string "default", i.e.
	Config("default") which sets the model to default settings.
	
	Keywords
	========
	Model : str
		Name of the internal field model to use, current models 
		available include:
		"jrm33" (default)|"jrm09"|"vip4"|"vit4"|"vipal"|"isaac"|
		"gsfc13ev"|"gsfc15ev"|"gsfc15evs"|"jpl15ev"|"jpl15evs"|
		"o4"|"o6"|"p11a"|"sha"|"u17ev"|"v117ev"|"none"	
	CartesianIn : bool
		If True (default) the inputs to the model will be expected to be 
		in Cartesian right-handed System III coordinates. If False, then
		the inputs should be in spherical polar coordinates.
	CartesianOut : bool
		If True (default) the output magnetic field will be in Cartesian
		right-handed System III coordinates. Otherwise, the magnetic 
		field components produced will be radial, meridional and 
		azimuthal.		
	Degree : int
		Maximum degree to use on the current model.
	'''

	#list the default arguments here
	defargs = {	'Model'			: 'jrm09',
				'CartesianIn'	: True,
				'CartesianOut'	: True,
				'Degree'		: 0}
				
	
	if len(args) == 1:
		if args[0] == "default":
			#set the configuration of the model to default
			_SetCFG(defargs)
			#note that we can still modify the default configuration
			#using the keywords provided after doing this		
	
	#return the current configuration
	cfg = _GetCFG()
		
						
	#check input kwargs
	#for those which exist  add
	#them to this object using the short name as the object tag
	#Otherwise use the existing value
		
	#the input keys
	ikeys = list(kwargs.keys())
		
	#current config keys
	ckeys = list(cfg.keys())
		
	#short and long name keys
	dkeys = list(defargs.keys())

	#check if the model has changed
	ChangeModel = False
	if 'Model' in ikeys and 'Model' in ikeys:
		if cfg['Model'] != kwargs['Model']:
			ChangeModel=True

	#loop through each one		
	for k in ckeys:
		if k in ikeys:
			#short name found in kwargs - add to this object
			cfg[k] = kwargs[k]
		else:
			#key not found, use existing parameter or default if model has changed
			if ChangeModel:
				cfg[k] = defargs[k]
			pass
		
	#check for additional keys and issue a warning
	for k in ikeys:
		if not ((k in dkeys) or (k in ckeys)):
			print("Keyword argument {:s} unrecognized, ignoring.".format(k))

	#update the configuration
	_SetCFG(cfg)
	
	#get a copy of the config stored in C++
	cfg = _GetCFG()

	
	return cfg
