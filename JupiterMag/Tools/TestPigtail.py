import numpy as np
import matplotlib.pyplot as plt



def TestPigtail(IntModel='jrm09',ExtModel='Con2020',
				t0=(302,12,50),t1=(302,20,30)):
	from ..TraceField import TraceField
	from ..Con2020 import Config
	from ..Internal._ReadTestPos import _ReadTestPos
	import DateTimeTools as TT
	

	#read the data
	print('Reading Data')
	data = _ReadTestPos()

	#get the time
	year = data.Year
	dayno = data.Doy
	

	
	#time limits (in dayno)
	d0 = t0[0] + t0[1]/24 + t0[2]/1440
	d1 = t1[0] + t1[1]/24 + t1[2]/1440
	
	
	#limit the dayno
	use = np.where((dayno >= d0) & (dayno <= d1))[0]
	data = data[use] 
	year = year[use]
	dayno = dayno[use]
	yr = data.Year[0]
	x = data.x
	y = data.y
	z = data.z

	#convert time
	dn = np.int32(dayno)
	ut = (dayno-dn)*24.0
	Date = TT.DayNotoDate(year,dn)

	
	#get the trace
	cfg = Config()
	Config(equation_type='analytic')

	T = TraceField(x,y,z,Verbose=True,IntModel=IntModel,ExtModel=ExtModel,Date=Date,ut=ut)

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

	s0 = '{:04d}-{:03d} {:02d}:{:02d}'.format(yr,t0[0],t0[1],t0[2])
	s1 = '{:04d}-{:03d} {:02d}:{:02d}'.format(yr,t1[0],t1[1],t1[2])

	ax.set_title(s0 + ' - ' + s1)

	return ax,T
