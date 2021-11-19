import numpy as np
from ._CFunctions import _CGetVIP4Config,_CSetVIP4Config
import ctypes as ct

def _GetCFG():
	'''
	Get the current config dictionary
	
	
	'''

	CartIn = np.zeros(1,dtype='bool')
	CartOut = np.zeros(1,dtype='bool')

	_CGetVIP4Config(CartIn,CartOut)
	
	cfg = {}
	cfg['CartesianIn'] = CartIn[0]
	cfg['CartesianOut'] = CartOut[0]
	return cfg

def _SetCFG(cfg):
	
	CartIn = np.array(cfg['CartesianOut'],dtype='bool')
	CartOut = np.array(cfg['CartesianOut'],dtype='bool')
	
	_CSetVIP4Config(CartIn,CartOut)
						
						
def Config(*args,**kwargs):
	'''
	Set and return the VIP4 model configuration
	
	'''

	#list the default arguments here
	defargs = {	'CartesianIn'	: True,
				'CartesianOut'	: True}
				
	
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

		
	#loop through each one		
	for k in ckeys:
		if k in ikeys:
			#short name found in kwargs - add to this object
			cfg[k] = kwargs[k]
		else:
			#key not found, use existing parameter
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
