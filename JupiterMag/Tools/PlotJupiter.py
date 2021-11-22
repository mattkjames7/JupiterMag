import numpy as np


def PlotJupiterXY(ax):
	'''
	Plot a circle
	
	'''
	
	
	#plot the outside edge
	a = np.arange(361)*np.pi/180.0
	x = np.cos(a)
	y = np.sin(a)
	ax.plot(x,y,color='black',linestyle='-',linewidth=2.0)
	

			
			
			
def PlotJupiterXZ(ax):
	'''
	Plot Jupiter in the X-Z (or Y or rho -Z plane, whatever)
	
	'''
	
	#Jupiter's equatorial and polar Radii (a and b, respectively)
	#in Rj, assuming Rj == equatorial radius
	t = np.arange(361.0)*np.pi/180.0
	
	a = 1.0
	b = 0.935

	x = a*np.cos(t)
	z = b*np.sin(t)
	
	ax.plot(x,z,color='black',linestyle='-',linewidth=2.0)
