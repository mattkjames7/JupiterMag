import numpy as np



def MagtoSIII(x,y,z,xt,xp):
	'''
	Convert from a dipole/current sheet based coordinate system back to
	Right-handed System III.
	
	Inputs
	======
	x : float
		x coordinate in dipole/current sheet coordinates
	y : float
		y coordinate in dipole/current sheet coordinates
	z : float
		z coordinate in dipole/current sheet coordinates
	xt : float
		Current sheet/dipole tilt (degrees)
	xp : float
		Azimuth of tilt (degrees)
		
	Returns
	=======
	ox : float
		x-coordinate 
	oy : float
		y-coordinate 
	oz : float
		z-coordinate 
	
	'''

	#some sines and cosines
	dtor = np.pi/180.0
	xtr = dtor*xt
	xpr = dtor*xp
	cosxt = np.cos(xtr)
	sinxt = np.sin(xtr)
	cosxp = np.cos(xpr)
	sinxp = np.sin(xpr)
	
	
	
