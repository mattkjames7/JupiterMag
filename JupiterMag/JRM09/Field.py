import numpy as np
from ._CFunctions import _CJRM09FieldArray,_CGetJRM09Config,_CSetJRM09Config

def Field(p0,p1,p2,CartIn=True,CartOut=True):
	'''
	Return the Con 2020 magnetic field vector(s). Check the model 
	config using JupiterMag.Con2020.Config() to see whether Cartesian or
	polar coordinates are used for input/output.
	
	Inputs
	======
	p0 : float
		Array/scalar containing x or r right-handed System III 
		coordinate (depending on PolIn)
	p1 : float
		Array/scalar containing y or theta right-handed System III 
		coordinate
	p2 : float
		Array/scalar containing z or phi right-handed System III 
		coordinate

	Returns
	=======
	B0 : float
		Either Bx or Br in nT
	B1 : float
		Either By or Btheta in nT
	B2 : float
		Either Bz or Bphi in nT
	
	'''

	#make sure that the inputs are the correct type
	if (hasattr(p0,'__iter__') == False):
		_p0 = np.array([p0]).astype('float64')
	else:
		_p0 = np.array(p0).astype('float64')
	if (hasattr(p1,'__iter__') == False):
		_p1 = np.array([p1]).astype('float64')
	else:
		_p1 = np.array(p1).astype('float64')
	if (hasattr(p2,'__iter__') == False):
		_p2 = np.array([p2]).astype('float64')
	else:
		_p2 = np.array(p2).astype('float64')
	_l = np.int32(np.size(_p0))

	_B0 = np.zeros(_l,dtype='float64')
	_B1 = np.zeros(_l,dtype='float64')
	_B2 = np.zeros(_l,dtype='float64')
	

	#call the model
	_CJRM09FieldArray(_l,_p0,_p1,_p2,_B0,_B1,_B2)
	
	return _B0,_B1,_B2
	
	
