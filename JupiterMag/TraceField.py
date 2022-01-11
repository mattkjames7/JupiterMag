import numpy as np
from ._CFunctions import _CTraceField
import ctypes
from ._ptr2D import _ptr2D
from . import Con2020
import matplotlib.pyplot as plt
from .Tools.PlotJupiter import PlotJupiterXY,PlotJupiterXZ
from .Tools.GetLegendHandLab import GetLegendHandLab
from . import Internal		

class TraceField(object):
	'''
	Object which stores the result of a magnetic field trace or a series 
	of traces performed using a combination of internal and external field
	models
	
	
	'''
	
	def __init__(self,x0,y0,z0,IntModel='jrm09',ExtModel='Con2020', 
				FlattenSingleTraces=True,**kwargs):
		'''
		Traces along the magnetic field given a starting set of 
		coordinates (or for multiple traces, arrays of starting 
		coordinates).
		
		Inputs
		=======
		x0: float
			scalar or array containing the x component of the starting 
			point(s).
		y0 : float
			scalar or array containing the y component of the starting 
			point(s).
		z0 : float
			scalar or array containing the z component of the starting 
			point(s).

		FlattenSingleTraces	: bool
			When set to True and if performing only a single trace to 
			flatten all of the arrays (position, magnetic field, etc.)
		Verbose	: bool
			Boolean, if True will display an indication of the progress 
			made during traces.
		TraceDir : int|str
			if set to 0 or 'both' then the trace will run in both 
			directions. Set to 1 to trace along the field direction
			(from south to north), or set to -1 to trace in the opposite
			direction to the magnetic field (north to south).
		
		Keyword arguments
		=================

		
		Model Fields
		============

				
		'''
		
		

		#Convert input variables to appropriate numpy dtype:
		self.x0 = np.array([x0]).flatten().astype("float64")
		self.y0 = np.array([y0]).flatten().astype("float64")
		self.z0 = np.array([z0]).flatten().astype("float64")
		self.n = np.int32(self.x0.size)

		self.IntModel = IntModel
		self.IntModelCode = ctypes.c_char_p(IntModel.encode('utf-8'))
		self.ExtModel = ExtModel
		self.ExtModelCode = ctypes.c_char_p(ExtModel.encode('utf-8'))

		#make sure models are in Cartesian
		Models = [IntModel.lower(),ExtModel]

		
		if ExtModel == "Con2020":
			Con2020.Config(CartesianIn=True,CartesianOut=True)
			

		#kwargs
		defargs = {	'MaxLen'	:		1000,
					'MaxStep'	:		1.0,
					'InitStep'	:		0.1,
					'MinStep'	:		0.0001,
					'ErrMax'	:		0.0001,
					'Delta'		:		0.05,
					'Verbose'	:		False,
					'TraceDir'	:		'both',
					'alpha' 	:		[]}
		dkeys = list(defargs.keys())
		kkeys = list(kwargs.keys())
		cfg = {}
		for k in dkeys:
			if k in kkeys:
				cfg[k] = kwargs[k]
			else:
				cfg[k] = defargs[k]
				
		
		self.Verbose = np.bool(cfg['Verbose'])
		self.MaxLen = np.int32(cfg['MaxLen'])
		self.MaxStep = np.float64(cfg['MaxStep'])
		self.InitStep = np.float64(cfg['InitStep'])
		self.MinStep = np.float64(cfg['MinStep'])
		self.ErrMax = np.float64(cfg['ErrMax'])
		self.Delta = np.float64(cfg['Delta'])
		TraceDir = cfg['TraceDir']
		if TraceDir == 'both':
			TraceDir = 0
		self.TraceDir = np.int32(TraceDir)



		self.x = np.zeros((self.n,self.MaxLen),dtype="float64") + np.nan
		self.y = np.zeros((self.n,self.MaxLen),dtype="float64") + np.nan
		self.z = np.zeros((self.n,self.MaxLen),dtype="float64") + np.nan
		self.Bx = np.zeros((self.n,self.MaxLen),dtype="float64") + np.nan
		self.By = np.zeros((self.n,self.MaxLen),dtype="float64") + np.nan
		self.Bz = np.zeros((self.n,self.MaxLen),dtype="float64") + np.nan

		self.nstep = np.zeros(self.n,dtype="int32")

		self.s = np.zeros((self.n,self.MaxLen),dtype="float64") + np.nan
		self.R = np.zeros((self.n,self.MaxLen),dtype="float64") + np.nan
		self.Rnorm = np.zeros((self.n,self.MaxLen),dtype="float64") + np.nan

		alpha = cfg['alpha']
		self.nalpha = np.int32(np.size(alpha))
		self.alpha = np.array(alpha).astype('float64')
		self.halpha = np.zeros((self.n*self.MaxLen*self.nalpha,),dtype="float64") + np.nan #hopefully this will be reshaped to (n,nalpha,MaxLen)
		self.FP = np.zeros((self.n,7),dtype="float64")

		_x = _ptr2D(self.x)
		_y = _ptr2D(self.y)
		_z = _ptr2D(self.z)

		_Bx = _ptr2D(self.Bx)
		_By = _ptr2D(self.By)
		_Bz = _ptr2D(self.Bz)
	
		
		_s = _ptr2D(self.s)
		_R = _ptr2D(self.R)
		_Rnorm = _ptr2D(self.Rnorm)		
		_FP = _ptr2D(self.FP)
		
		#call the C code
		_CTraceField(	self.n,self.x0,self.y0,self.z0,
						self.IntModelCode,self.ExtModelCode,
						self.MaxLen,self.MaxStep,self.InitStep,
						self.MinStep,self.ErrMax,self.Delta,
						self.Verbose,self.TraceDir,
						self.nstep,
						_x,_y,_z,
						_Bx,_By,_Bz,
						_R,_s,_Rnorm,_FP,
						self.nalpha,self.alpha,self.halpha)

		#reshape the footprints
		fpnames = ['LatN','LonN','LatS','LonS','LonEq','Rmax','FlLen']


		
		#flatten things and unpack footprints
		if self.n == 1 and FlattenSingleTraces:
			flat = ['nstep','x','y','z','Bx','By','Bz','s','R','Rnorm']
			for f in flat:
				self.__dict__[f] = self.__dict__[f][0]
			self.halpha = (self.halpha.reshape((self.n,self.nalpha,self.MaxLen)))[0]
			for i in range(0,7):
				setattr(self,fpnames[i],self.FP[0,i])
		else:
			self.halpha = self.halpha.reshape((self.n,self.nalpha,self.MaxLen))
			for i in range(0,7):
				setattr(self,fpnames[i],self.FP[:,i])

	def PlotXZ(self,ind='all',fig=None,maps=[1,1,0,0],label=None,color='black'):
		'''
		plot field lines in the X-Z plane
		
		'''
		
		if ind == 'all':
			ind = np.arange(self.n)
		elif np.size(ind) == 1:
			ind = np.array([ind]).flatten()
		else:
			ind = np.array(ind)
			
		
		if fig is None:
			fig = plt
			fig.figure()
		if hasattr(fig,'Axes'):	
			ax = fig.subplot2grid((maps[1],maps[0]),(maps[3],maps[2]))
		else:
			ax = fig
		
		if np.size(np.shape(self.x)) == 2:
			x = self.x[ind].T
			z = self.z[ind].T
		else:
			x = self.x
			z = self.z
			

		ln = ax.plot(x,z,color=color)
		if not label is None:
			hs,ls = GetLegendHandLab(ax)
			hs.append(ln[0])
			ls.append(label)
			ax.legend(hs,ls)
		
		ax.set_ylabel('$z_{SIII}$ (R$_J$)')
		ax.set_xlabel('$x_{SIII}$ (R$_J$)')

		mxx = np.nanmax(x)
		mxz = np.nanmax(z)
		mx = 1.1*np.nanmax([mxx,mxz])		
		ax.set_xlim(-mx,mx)
		ax.set_ylim(-mx,mx)
		
		PlotJupiterXZ(ax)
		ax.set_aspect(1.0)

		return ax
	
	def PlotXY(self,ind='all',fig=None,maps=[1,1,0,0],label=None,color='black'):
		'''
		plot field lines in the X-Y plane
		
		'''
		
		if ind == 'all':
			ind = np.arange(self.n)
		elif np.size(ind) == 1:
			ind = np.array([ind]).flatten()
		else:
			ind = np.array(ind)
			
		
		if fig is None:
			fig = plt
			fig.figure()
		if hasattr(fig,'Axes'):	
			ax = fig.subplot2grid((maps[1],maps[0]),(maps[3],maps[2]))
		else:
			ax = fig
		
		if np.size(np.shape(self.x)) == 2:
			x = self.x[ind].T
			y = self.y[ind].T
		else:
			x = self.x
			y = self.y
			
		ln = ax.plot(y,x,color=color)
		if not label is None:
			hs,ls = GetLegendHandLab(ax)
			hs.append(ln[0])
			ls.append(label)
			ax.legend(hs,ls)
		yl = ax.get_xlim()
		ax.set_xlim(yl[::-1])
		
		ax.set_xlabel('$y_{SIII}$ (R$_J$)')
		ax.set_ylabel('$x_{SIII}$ (R$_J$)')

		mxx = np.nanmax(x)
		mxy = np.nanmax(y)
		mx = 1.1*np.nanmax([mxx,mxy])		
		ax.set_xlim(mx,-mx)
		ax.set_ylim(-mx,mx)
		
		PlotJupiterXY(ax)
		ax.set_aspect(1.0)
		return ax
	
	def PlotRhoZ(self,ind='all',fig=None,maps=[1,1,0,0],label=None,color='black'):
		'''
		plot field lines in the rho-Z plane
		
		'''
		
		if ind == 'all':
			ind = np.arange(self.n)
		elif np.size(ind) == 1:
			ind = np.array([ind]).flatten()
		else:
			ind = np.array(ind)
			
		
		if fig is None:
			fig = plt
			fig.figure()
		if hasattr(fig,'Axes'):	
			ax = fig.subplot2grid((maps[1],maps[0]),(maps[3],maps[2]))
		else:
			ax = fig
		
		if np.size(np.shape(self.x)) == 2:
			x = self.x[ind].T
			y = self.y[ind].T
			z = self.z[ind].T
		else:
			x = self.x
			y = self.y
			z = self.z
		
		r = np.sqrt(x**2 + y**2)
		ln = ax.plot(r,z,color=color)
		if not label is None:
			hs,ls = GetLegendHandLab(ax)
			hs.append(ln[0])
			ls.append(label)
			ax.legend(hs,ls)
		
		ax.set_ylabel('$z_{SIII}$ (R$_J$)')
		ax.set_xlabel(r'$\rho_{SIII}$ (R$_J$)')

		mxr = np.nanmax(r)
		mxz = np.nanmax(z)
		mx = 1.1*np.nanmax([mxr,mxz])		
		ax.set_xlim(-mx,mx)
		ax.set_ylim(-mx,mx)
		
		PlotJupiterXZ(ax)
		ax.set_aspect(1.0)
		return ax
	
	
	def PlotHalpha(self,TI='all',AI='all',fig=None,maps=[1,1,0,0]):
		'''
		Plot h_alpha (see Singer et al 1982) for a field line.
		
		Inputs
		======
		TI : int|str
			Index of trace to plot. TI='all' will plot for all traces.
		AI : int|str
			Index of alpha angle to plot for. AI will plot all alphas.
		fig : None|matplotlib.pyplot|matplotlib.pyplot.Axes
			None - a new figure will be created with new axes
			matplotlib.pyplot - existing figure, new axes
			matplotlib.pyplot.Axes - existing axes instance to be used
				(maps ignored in the case).
		maps : list|tuple|numpy.ndarray
			Four element array-like, denoting subplot position,
			e.g. [xmaps,ymaps,xmap,ymap]
				xmaps : number of subplots in x-direction
				ymaps : number of subplots in y-direction
				xmap : position index (0 is left)
				ymap : position index (0 is top)
		
		
		'''
		if AI == 'all':
			AI = np.arange(self.nalpha)
		
		if np.size(AI) == 1:
			AI = np.array([AI]).flatten()
			
		if TI == 'all':
			TI = np.arange(self.n)
		
		if np.size(TI) == 1:
			TI = np.array([TI]).flatten()
			
		if fig is None:
			fig = plt
			fig.figure()
		if hasattr(fig,'Axes'):	
			ax = fig.subplot2grid((maps[1],maps[0]),(maps[3],maps[2]))
		else:
			ax = fig		
			
		for t in TI:
			if np.size(np.shape(self.x)) == 1:
				#trace arrays have been flattened 
				for a in AI:
					ax.plot(self.s,self.halpha[a],label=r'$\alpha=${:5.1f}'.format(self.alpha[a]))
			else:
				for a in AI:
					ax.plot(self.s[t],self.halpha[t,a],label=r'Trace {:d} $\alpha=${:5.1f}'.format(t,self.alpha[a]))

		ax.legend()
		ax.set_xlabel(r'$s$ (R$_J$)')
		ax.set_ylabel(r'$h_{\alpha}$')

		return ax
		
	def PlotPigtail(self,Proj='normal',ShowLabels=True,Time=None,
					Hemisphere='both',colatlim=None,
					fig=None,maps=[1,1,0,0],**kwargs):
		'''
		Pigtail plot.
		
		Inputs
		======
		Proj : str
			'normal' : plot footprints on latitude/longitude plot
			'abnormal' : plot as though we are looking down on the pole
		ShowLabels : bool
			This will display some sort of time axis, if Time is provided
		Time : None|float64|(int32,float32)
			Time of each trace - must have same number of elements as 
			there are traces.
			float64 : continuous time
			(int32,float32) : (Date formatted yyyymmdd,UT in hours)
		Hemisphere : str
			'north'|'south'|'both'

		
		'''
		
		#get the stuff to plot
		rn = np.abs(self.LatN)
		rs = np.abs(self.LatS)
		tn = self.LonN*np.pi/180.0
		ts = self.LonS*np.pi/180.0
		if Proj == 'abnormal':
			rn = np.sin(rn)
			rs = np.sin(rs)

		#lower latitude limit
		if colatlim is None:
			if Proj == 'normal':
				colatlim = np.min([rn.min(),rs.min()])
			else:
				colatlim = 1.0
		if Proj == 'normal':
			rlim = [90.0,colatlim]
		else:
			rlim = [0.0,1.0]
		
		if fig is None:
			fig = plt
			fig.figure()
		if hasattr(fig,'Axes'):	
			ax = fig.subplot2grid((maps[1],maps[0]),(maps[3],maps[2]),projection='polar')
		else:
			ax = fig		
		ax.set_theta_zero_location("N")
		ax.set_rlim(rlim)
		if Hemisphere.lower() in ['both','north']:
			ax.plot(tn,rn,linewidth=kwargs.get('linewidth',2.0),color=kwargs.get('color','red'),label='North')	
		if Hemisphere.lower() in ['both','south']:
			ax.plot(ts,rs,linewidth=kwargs.get('linewidth',2.0),color=kwargs.get('color','orange'),label='South')	
		
		ax.legend()
		
		return ax
		
