import numpy as np
import matplotlib.pyplot as plt


def TestTrace(IntModel='jrm33',ExtModel='Con2020',fig=None,maps=[1,1,0,0],color='green'):
	from ..TraceField import TraceField
	from ..Con2020 import Config
	
	#set the starting coords
	n = 7
	x = np.linspace(2.0,30.0,n)
	x = np.append(-x[::-1],x)
	y = np.zeros(n*2)
	z = np.zeros(n*2)
	
	#get the trace
	cfg = Config()
	Config(equation_type='analytic')
	T = TraceField(x,y,z,Verbose=True,IntModel=IntModel,ExtModel=ExtModel)
	Config(cfg)

	#plot it
	lab = ''
	if not IntModel.upper() == 'NONE':
		lab += IntModel
	
	if not ExtModel.upper() == 'NONE':
		if not lab == '':
			lab += ' + '
		lab += ExtModel
	
	ax = T.PlotXZ(fig=fig,maps=maps,label=lab,color=color)


	return ax


def CompareTrace():
	from ..TraceField import TraceField
	from ..Con2020 import Config
		

	#get some starting coords
	n = 8
	theta = (180.0 - np.linspace(22.5,35,n))*np.pi/180.0
	r = np.ones(n)
	x = r*np.sin(theta)
	y = np.zeros(n)
	z = r*np.cos(theta)
	
	#get traces with and without the external field
	cfg = Config()
	Config(equation_type='analytic')	
	T0 = TraceField(x,y,z,Verbose=True,IntModel='jrm33',ExtModel='none')
	T1 = TraceField(x,y,z,Verbose=True,IntModel='jrm33',ExtModel='Con2020')
	Config(cfg)
	
	#plot them
	ax = T0.PlotRhoZ(label='JRM33',color='black')
	ax = T1.PlotRhoZ(fig=ax,label='JRM33 + Con2020',color='red')
	
	ax.set_xlim(-2.0,15.0)
	ax.set_ylim(-6.0,6.0)
	
	return ax
