import numpy as np
from .Field import Field
from .Config import Config
import time
from ..Con2020._ReadTestData import _ReadTestData
import matplotlib.pyplot as plt

def Test(Model='JRM09',R=0.85,MaxDeg=None,Comp='r',scale = [-60.0,60.0],
			levels=np.linspace(-50,50,11),fig=None,maps=[1,1,0,0]):
	'''
	This is a simple function to test the model by recreating a plot in 
	Connerney et al 2018 (figure 4, sort of).
	
	Inputs
	======
	R : float
		The radial distance to evaluate the model at.
	MaxDeg : int
		Maximum model degree to calculate.
	
	'''
	try:
		import matplotlib.pyplot as plt
		import matplotlib.colors as colors
		from mpl_toolkits.axes_grid1 import make_axes_locatable
	except:
		raise SystemError('This function requires "matplotlib" to be instaled')

	#get the coordinates to calculate the model at
	lat = np.linspace(-90,90,181)
	lon = np.linspace(0.0,360.0,361)
	latc = 0.5*(lat[1:] + lat[:-1])
	lonc = 0.5*(lon[1:] + lon[:-1])
	long,latg = np.meshgrid(lon,lat)
	longc,latgc = np.meshgrid(lonc,latc)

	longcr = longc*np.pi/180.0
	latgcr = (90.0 - latgc)*np.pi/180.0
	r = np.zeros(longcr.shape) + R
	
	#get original config
	cfg0 = Config()
	
	#set to polar
	Config(Model=Model,CartesianIn=False,CartesianOut=False)
	
	#calculate the model
	Br,Bt,Bp = Field(r,latgcr,longcr,MaxDeg=MaxDeg)
	
	#B = np.sqrt(Br**2 + Bt**2 + Bp**2)
	if Comp == 'r':
		B = Br
		Blab = '$B_r$'
	elif Comp == 'p':
		B = Bp
		Blab = '$B_p$'
	elif Comp == 't':
		B = Bt
		Blab = '$B_t$'
	else:
		B = np.sqrt(Br**2 + Bt**2 + Bp**2)
		Blab = '$|B|$'
	
	#convert to Gauss
	Bg = B.reshape(longcr.shape)*1e-5
	
	#restore config
	Config(**cfg0)

	
	#set norm
	norm = colors.Normalize(vmin=scale[0],vmax=scale[1])	
		
	if fig is None:
		fig = plt
		fig.figure()
	if hasattr(fig,'Axes'):	
		ax = fig.subplot2grid((maps[1],maps[0]),(maps[3],maps[2]))
	else:
		ax = fig

	ax.set_aspect(1.0)
	ax.set_xlabel('SIII East Longitude ($^\circ$)')
	ax.set_ylabel('SIII Latitude ($^\circ$)')
	if MaxDeg is None:
		ax.set_title(Model)
	else:
		ax.set_title(Model + ' (Deg={:d})'.format(MaxDeg))
	sm = ax.pcolormesh(long,latg,Bg,cmap='RdYlBu_r',norm=norm)
	ct = ax.contour(longc,latgc,Bg,colors='grey',levels=levels)
	ax.clabel(ct, inline=True, fontsize=8,fmt='%2d')

	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="5%", pad=0.05)

	cbar = plt.colorbar(sm,cax=cax) 
	cbar.set_label('{:s} (Gauss) at $r$ = {:4.2f}'.format(Blab,R) + ' R$_{J}$')
	
	return ax

def TimeVector():

	#get original config
	cfg0 = Config()
	
	#set to polar
	Config(CartesianIn=False,CartesianOut=False)
	
	#read data
	data = _ReadTestData()
	
	#call model
	t0 = time.time()
	Br,Bt,Bp = Field(data.r,data.t,data.p)
	t1 = time.time()
	
	#restore config
	Config(**cfg0)
	
	print('Time: {:f}s'.format(t1-t0))
	
def TimeScalar():

	#get original config
	cfg0 = Config()
	
	#set to polar
	Config(CartesianIn=False,CartesianOut=False)
	
	#read data
	data = _ReadTestData()
	
	#call model
	t0 = time.time()
	for i in range(0,data.size):
		Br,Bt,Bp = Field(data.r[i],data.t[i],data.p[i])
	t1 = time.time()
	
	#restore config
	Config(**cfg0)
	
	print('Time: {:f}s'.format(t1-t0))


def TestOutput(fname=None,Model='VIP4',MaxDeg=None):

	#get original config
	cfg0 = Config()
	
	#set to polar
	Config(Model=Model,CartesianIn=False,CartesianOut=False)
	
	#positions to test
	r = np.array([3,3,3,3, 3,3,3,3, 3,3,3,3, 3,3,3,3],dtype='float64')
	theta = np.array([10,10,10,10,55,55,55,55,90,90,90,90,130,130,130,130],dtype='float64')
	phi = np.array([0,27,180,340, 0,27,180,340, 0,27,180,340, 0,27,180,340],dtype='float64')
	
	#model output
	Br,Bt,Bp = Field(r,theta*np.pi/180.0,phi*np.pi/180.0,MaxDeg=MaxDeg)

	
	#restore config
	Config(**cfg0)

	#save to file
	lines = []
	
	out = '  R  | Theta |  Phi  |         Br         |         Bt         |         Bp         '
	print(out)
	lines.append(out)

	out = '-----|-------|-------|--------------------|--------------------|--------------------' 
	print(out)
	lines.append(out)
	for i in range(0,r.size):
		out = ' {:3.1f} | {:5.1f} | {:5.1f} | {:18.11f} | {:18.11f} | {:18.11f}'.format(r[i],theta[i],phi[i],Br[i],Bt[i],Bp[i])
		print(out)
		lines.append(out)
	
	if not fname is None:
		f = open(fname,'w')
		for l in lines:
			f.write(l+'\n')
		f.close()


def JRM33Fig5(MaxDeg=13,fig=None,maps=[1,1,0,0]):
	Test('JRM33',MaxDeg=MaxDeg,Comp='B',scale=[0,20],levels=np.linspace(0,20,11),R=1.0,fig=fig,maps=maps)
	plt.savefig('JRM33-Fig5-Deg{:d}.png'.format(MaxDeg))

def JRM33Fig7(MaxDeg=13,fig=None,maps=[1,1,0,0]):
	Test('JRM33',MaxDeg=MaxDeg,Comp='r',scale=[-80,80],levels=np.linspace(-80,80,17),R=0.85,fig=fig,maps=maps)
	plt.savefig('JRM33-Fig7-Deg{:d}.png'.format(MaxDeg))
