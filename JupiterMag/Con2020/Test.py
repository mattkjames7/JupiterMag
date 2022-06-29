import numpy as np
import matplotlib.pyplot as plt
from .Field import Field
from .Config import Config
from ._ReadTestData import _ReadTestData

def _RCrossings(doy,r,r0,r1):
	'''
	Return the day numbers where r0 and r1 were crossed.
	
	Inputs
	======
	doy : float
		Day number array.
	r : float
		Array of R.
	r0 : float
		Inside edge of the current.
	r0 : float 
		Outside edge of the current.
		
	Returns
	=======
	dc : float
		Day numbers of the r0/r1 crossings.
	
	'''
	
	rc0 = ((r[1:] >= r0) & (r[:-1] < r0)) | ((r[:-1] >= r0) & (r[1:] < r0))
	rc1 = ((r[1:] >= r1) & (r[:-1] < r1)) | ((r[:-1] >= r1) & (r[1:] < r1))
	rc0 = np.where(rc0)[0]
	rc1 = np.where(rc1)[0]
	rc = np.append(rc0,rc1)
	rc.sort()
	
	dc = 0.5*(doy[rc] + doy[rc+1])
	
	return dc
	
	
	
	
def Test():
	'''
	Run a quick test to see if the model works.
	
	'''

	#read the data
	print('Reading Data')
	data = _ReadTestData()

	#get the time
	year = data.Year
	dayno = data.Day
	
	#limit the dayno
	use = np.where((dayno >= 290) & (dayno <= 315))[0]
	data = data[use] 


	#call the model code
	print('Calling Model')
	cfg = Config('default',CartesianIn=False,CartesianOut=False,Edwards=False)
	Br,Bt,Bp = Field(data.r,data.t,data.p)
	
	#get the r0/r1 crossings
	dc = _RCrossings(data.Day,data.r,cfg['r0'],cfg['r1'])

	#create a plot window
	plt.figure(figsize=(11,8))
	
	#create the subplots
	ax0 = plt.subplot2grid((4,1),(0,0))
	ax1 = plt.subplot2grid((4,1),(1,0))
	ax2 = plt.subplot2grid((4,1),(2,0))
	ax3 = plt.subplot2grid((4,1),(3,0))
	
	#plot each component
	ax0.plot(data.Day,Br,color='k',label=r'$B_{r}$ (nT)')
	ax1.plot(data.Day,Bt,color='k',label=r'$B_{\theta}$ (nT)')
	ax2.plot(data.Day,Bp,color='k',label=r'$B_{\phi}$')
	ax3.plot(data.Day,data.r,color='k',label=r'$r$')
	
	#fix y limits
	y0 = ax0.get_ylim()
	y1 = ax1.get_ylim()
	y2 = ax2.get_ylim()
	y3 = ax3.get_ylim()
	ax0.set_ylim(y0)
	ax1.set_ylim(y1)
	ax2.set_ylim(y2)
	ax3.set_ylim(y3)
	
	#and x limits
	ax0.set_xlim([290,315])
	ax1.set_xlim([290,315])
	ax2.set_xlim([290,315])
	ax3.set_xlim([290,315])
	
	#plot r0/r1 crossings
	ax0.vlines(dc,y0[0],y0[1],color='k',linestyle='--')
	ax1.vlines(dc,y1[0],y1[1],color='k',linestyle='--')
	ax2.vlines(dc,y2[0],y2[1],color='k',linestyle='--')
	ax3.vlines(dc,y3[0],y3[1],color='k',linestyle='--')
	
	#y labels plot labels
	ax0.set_ylabel(r'$B_r$ / nT')
	ax1.set_ylabel(r'$B_{\theta}$ / nT')
	ax2.set_ylabel(r'$B_{\phi}$ / nT')
	ax3.set_ylabel(r'$r$ / R$_J$')

	#title
	ax0.set_title(r'PJ16, con2020: $\mu_0I_{MD}/2$=' + '{:5.1f}, $R_0$={:3.1f} R$_J$, $R_1$={:4.1f} R$_J$'.format(cfg['i_rho'],cfg['r0'],cfg['r1']))
	
	#x labels
	ax0.set_xticks([])
	ax1.set_xticks([])
	ax2.set_xticks([])
	ax3.set_xlabel('DOY (2018)')
	
	plt.subplots_adjust(hspace=0.0)
