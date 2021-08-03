import numpy as np
from ._CFunctions import _CCon2020Field

def Con2020Field(p0,p1,p2,PolIn=False,PolOut=False):
	'''
	Call one of the Con2020 field model for Jupiter.
	
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
	PolIn : bool
		If True - input coordinates are r, theta and phi (radial 
		distance in Rj, colatitude in radians and east longitude in
		radians), otherwise inputs are Cartesian x, y and z in Rj.
	PolOut : bool
		If True - output field vector is in polar coordinates, otherwise
		it is returned in Cartesian System III coordinates.

	Returns
	=======
	B0 : float
		Either Bx or Br (depending on PolOut) in nT
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
	_PolIn = np.bool8(PolIn)
	_PolOut = np.bool8(PolOut)
	_B0 = np.zeros(_l,dtype='float64')
	_B1 = np.zeros(_l,dtype='float64')
	_B2 = np.zeros(_l,dtype='float64')
	
	#call the model
	_CCon2020Field(_l,_p0,_p1,_p2,_B0,_B1,_B2,_PolIn,_PolOut)
	
	return _B0,_B1,_B2
