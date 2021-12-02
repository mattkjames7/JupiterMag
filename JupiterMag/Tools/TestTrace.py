import numpy as np
import matplotlib.pyplot as plt


def TestTrace(IntModel='JRM09',ExtModel='Con2020',fig=None,maps=[1,1,0,0],color='green'):
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
