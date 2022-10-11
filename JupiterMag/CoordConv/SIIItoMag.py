import numpy as np



def SIIItoMag(x,y,z,xt,xp):
	'''
	Convert Right-handed System III coordinates to a dipole/current
	sheet based coordinate system.
	
	Inputs
	======
	x : float
		x coordinate in SIII
	y : float
		y coordinate in SIII
	z : float
		z coordinate in SIII
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
	
	#now for rotating the coordinates
	xtmp = x*cosxp + y*sinxp
	ox = xtmp*cosxt + z*sinxt
	oy = y*cosxp - x*sinxp
	oz = z*cosxt - xtmp*sinxt
	
	return ox,oy,oz
	
	
