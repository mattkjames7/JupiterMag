import numpy as np
import matplotlib.pyplot as plt


def TestTrace():
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
	T = TraceField(x,y,z,Verbose=True)
	Config(cfg)

	#plot it
	ax = T.PlotXZ()


	return ax
