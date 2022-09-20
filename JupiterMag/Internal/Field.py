import numpy as np
from ._CFunctions import _CInternalField,_CInternalFieldDeg

def Field(p0,p1,p2,MaxDeg=None):
	'''
	Return the internal magnetic field vector(s). Check the model 
	config using JupiterMag.Internal.Config() to see whether Cartesian or
	polar coordinates are used for input/output and to set the model.
	
	Inputs
	======
	p0 : float
		Array/scalar containing x or r right-handed System III 
		coordinate 
	p1 : float
		Array/scalar containing y or theta right-handed System III 
		coordinate
	p2 : float
		Array/scalar containing z or phi right-handed System III 
		coordinate
	MaxDeg : None|int
		Maximum model degree to use. If None then the default value
		(model dependant) will be used.

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
	_p0 = np.array(p0,dtype='float64')
	_p1 = np.array(p1,dtype='float64')
	_p2 = np.array(p2,dtype='float64')
	_l = np.int32(np.size(_p0))

	_B0 = np.zeros(_l,dtype='float64')
	_B1 = np.zeros(_l,dtype='float64')
	_B2 = np.zeros(_l,dtype='float64')
	

	#call the model
	if MaxDeg is None:
		_CInternalField(_l,_p0,_p1,_p2,_B0,_B1,_B2)
	else:
		_MaxDeg = np.int32(MaxDeg)
		_CInternalFieldDeg(_l,_p0,_p1,_p2,_MaxDeg,_B0,_B1,_B2)
	
	return _B0,_B1,_B2
	
	
