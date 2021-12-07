import numpy as np
import matplotlib.pyplot as plt



def TestPigtail(IntModel='JRM09',ExtModel='Con2020'):
	from ..TraceField import TraceField
	from ..Con2020 import Config
	from ..Con2020._ReadTestData import _ReadTestData

	#read the data
	print('Reading Data')
	data = _ReadTestData()

	#get the time
	year = data.Year
	dayno = data.Day
	
	
	#limit the dayno
	use = np.where((dayno >= 302.534) & (dayno <= 302.8736))[0]
	data = data[use] 
	r = data.r
	t = data.t
	p = data.p
	rho = r*np.sin(t)
	z = r*np.cos(t)
	x = rho*np.cos(p)
	y = rho*np.sin(p)
	
	
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
	
	ax = T.PlotPigtail()


	return ax,T
