import numpy as np
from ._CFunctions import _CTraceField
import ctypes
from ._ptr2D import _ptr2D
from . import Con2020
import matplotlib.pyplot as plt
from .Tools.PlotJupiter import PlotJupiterXY,PlotJupiterXZ
from .Tools.GetLegendHandLab import GetLegendHandLab
from . import Internal		
import DateTimeTools as TT
from scipy.interpolate import interp1d
from .Tools.JupiterOval import JupiterOvalNorth,JupiterOvalSouth
import PyFileIO as pf

class TraceField(object):
	'''
	Object which stores the result of a magnetic field trace or a series 
	of traces performed using a combination of internal and external field
	models
	
	
	'''
	def __init__(self,*args,**kwargs):
		
		#check if we are loading from file, or creating new traces
		if len(args) == 1:
			#read from file or dict
			if isinstance(args[0],dict):
				#assume that the dictionary provided is a TraceField dict
				self.__dict__ = args[0]
			else:
				#load from file
				self._Load(*args)
		elif len(args) == 3:
			#new traces
			self._Trace(*args,**kwargs)
		else:
			#something's wrong
			print('TraceField was supplied with {:d} arguments...'.format(len(args)))
			print('Either use 1 string (file name), or')
			print('use 3 inputs (x,y,z)')
			return None
		
	
	
	
	def _Load(self,fname):
		self.__dict__ = pf.LoadObject(fname)
	
	
	
	def _Trace(self,x0,y0,z0,IntModel='jrm33',ExtModel='Con2020',**kwargs):
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
		InternalModel : str
			Name of the internal field model to use, current models 
			available include:
			"jrm33" (default)|"jrm09"|"vip4"|"vit4"|"vipal"|"isaac"|
			"gsfc13ev"|"gsfc15ev"|"gsfc15evs"|"jpl15ev"|"jpl15evs"|
			"o4"|"o6"|"p11a"|"sha"|"u17ev"|"v117ev"|"none"
		ExtModel : str
			External field model, currently only:
			"Con2020"|"none"		

		Keyword arguments
		=================
		Verbose	: bool
			Boolean, if True will display an indication of the progress 
			made during traces.
		TraceDir : int|str
			if set to 0 or 'both' then the trace will run in both 
			directions. Set to 1 to trace along the field direction
			(from south to north), or set to -1 to trace in the opposite
			direction to the magnetic field (north to south).
		MaxLen : int
			Maximum total number of trace steps
		MaxStep : float
			Length of maximum step size in planetary radii.
		InitStep : float
			Initial step size in planetary radii
		MinStep : float
			Minimum step size in planetary radii
		ErrMax : float
			Maximum allowed error in Runge-Kutta algorithm.
		alpha : float
			Array-like list of polarization angles for which to 
			calculate h_alpha (see Singer et al., 1981,
			doi: 10.1029/JA086iA06p04589)
		Delta : This is the separation between the equatorial footprints
			of the field lines used to calculate h_alpha.
			
		Member Functions
		================
		PlotXY()		Plots the field traces in the X-Y plane
		PlotXZ()		Plots the field traces in the X-Z plane
		PlotRhoZ()		Plots the field traces in the Rho-Z plane
		PlotHalpha()	Plots h_alpha along a field line
		PlotPigtail()	Fails to plot the pigtail plots.
		
		Member Variables
		================
		nstep 			Number of steps for each trace

		The following variables either have shape (n,MaxLen) or 
		(MaxLen,) if a single trace has been flattened. The elements of
		each trace past nstep[i] are filled with NANs.

		x				Trace x position 
		y				Trace y position
		z				Trace z position
		Bx				Trace field
		By				Trace field
		Bz				Trace field
		R				Radial distance along field line
		Rnorm			R/Rmax
		s				Distance along field line 
		
		These variables describe things such as footprints
		
		LatN			Latitude of northern footprints (degrees)
		LonN			Longitude of northern footprints (degrees)
		LatS			Latitude of southern foorprints (degrees)
		LonS			Longitude of southern footprints (degrees)
		LonEq			Longitude of magnetic equatorial footprint 
						(degrees)
		Rmax			Radial distance of the furthest point along the
						field line (planetary radii)
		FlLen			Length of field lines (planetary radii)
		
		
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
		self.nExt,self.ExtModelCode = self._WrapExtFuncs(ExtModel)

		#make sure models are in Cartesian
		Models = [IntModel.lower(),ExtModel]

		if ExtModel == "Con2020":
			Con2020.Config(CartesianIn=True,CartesianOut=True)
			
		#check if time has been supplied
		self.Time = False
		if 'Time' in kwargs:
			self._StoreTime(kwargs['Time'])
		elif 'Date' in kwargs and 'ut' in kwargs:
			self._StoreTime((kwargs['Date'],kwargs['ut']))


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
				
		
		self.Verbose = np.bool8(cfg['Verbose'])
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
		self.FP = np.zeros((self.n,49),dtype="float64")

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
						self.IntModelCode,self.nExt,self.ExtModelCode,
						self.MaxLen,self.MaxStep,self.InitStep,
						self.MinStep,self.ErrMax,self.Delta,
						self.Verbose,self.TraceDir,
						self.nstep,
						_x,_y,_z,
						_Bx,_By,_Bz,
						_R,_s,_Rnorm,_FP,
						self.nalpha,self.alpha,self.halpha)

		#reshape the footprints
		#fpnames = ['LatN','LonN','LatS','LonS','LonEq','Rmax','FlLen']


		
		#unpack footprints
		self._UnpackFootprints()

		#reshape h alpha
		self.halpha = self.halpha.reshape((self.n,self.nalpha,self.MaxLen))
		#for i in range(0,7):
			#setattr(self,fpnames[i],self.FP[:,i])


	def _UnpackFootprints(self):

		dtype = [	('xn3','float64'),
	    			('yn3','float64'),
					('zn3','float64'),
					('xs3','float64'),
	    			('ys3','float64'),
					('zs3','float64'),
					('xnm','float64'),
	    			('ynm','float64'),
					('znm','float64'),
					('xsm','float64'),
	    			('ysm','float64'),
					('zsm','float64'),
					('lonn','float64'),
					('latn','float64'),
					('mlonn','float64'),
					('mlatn','float64'),
					('lons','float64'),
					('lats','float64'),
					('mlons','float64'),
					('mlats','float64')]
		
		edtype = [	('x3','float64'),
	    			('y3','float64'),
					('z3','float64'),
					('xm','float64'),
	    			('ym','float64'),
					('zm','float64'),
					('lshell','float64'),
					('mlone','float64'),
					('fllen','float64')]

		iinds = np.array([0,1,2,3,4,5,6,7,8,9,10,11,30,31,32,33,34,35,36,37])
		sinds = np.array([12,13,14,15,16,17,18,19,20,21,22,23,38,39,40,41,42,43,44,45])
		einds = np.array([24,25,26,27,28,29,46,47,48])

		self.ionosphere = np.recarray(self.n,dtype=dtype)
		self.surface = np.recarray(self.n,dtype=dtype)
		self.equator = np.recarray(self.n,dtype=edtype)

		for i,f in enumerate(self.ionosphere.dtype.names):
			self.ionosphere[f] = self.FP[:,iinds[i]]


		for i,f in enumerate(self.surface.dtype.names):
			self.surface[f] = self.FP[:,sinds[i]]


		for i,f in enumerate(self.equator.dtype.names):
			self.equator[f] = self.FP[:,einds[i]]

	def _WrapExtFuncs(self,ExtFuncs):
		'''
		This will deal with the string/list of strings denoting the
		names of the external field functions to be used in the traces.
		It converts them into a compatible type with char**
		
		Inputs
		======
		ExtFuncs : str|list
			Name(s) of external field functions.
			
		Returns
		=======
		nExt : int32
			Number of external functions
		ExtPtr : ctypes.POINTER(ctypes.c_char_p)
		
		'''
		
		#convert to list if needed
		if isinstance(ExtFuncs,str):
			ExtFuncs = [ExtFuncs]
		
		#get the length
		nExt = np.int32(len(ExtFuncs))
		
		#create the pointer
		ExtPtr = (ctypes.c_char_p*nExt)()
		
		#encode the strings as bytes
		for i in range(0,nExt):
			ExtPtr[i] = ExtFuncs[i].encode('utf-8')
			
		return nExt,ExtPtr


	def PlotXZ(self,ind='all',fig=None,maps=[1,1,0,0],label=None,color='black'):
		'''
		Plot field lines in the X-Z plane
		
		Inputs
		======
		ind : int|str
			Index of trace to plot. Can be scalar or an array. If set 
			ind='all' then all traces will be plotted.
		fig : None|pyplot|pyplot.Axes instance
			None - new figure will be created
			pyplot - new subplot will be created on existing figure
			pyplot.Axes - existing subplot will be used
		maps : list
			4-element array-like to determine the subplot position,
			ignored when fig=pyplot.Axes.
			maps = [xmaps,ymaps,xmap,ymap]
			xmaps - number of subplots in x-direction
			ymaps - number of subplots in y-direction
			xmap - x position of this subplot
			ymap - y position of this subplot
		label : None|str
			Add label to traces.
		color : str|array-like
			Colour to plot the field lines
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
		
		x = self.x[ind]
		z = self.z[ind]
	
		mx = 1.5
		for i in range(0,ind.size):
			ln = ax.plot(x[i],z[i],color=color)
			mx = np.nanmax([mx,np.nanmax(np.abs(x[i])),np.nanmax(np.abs(z[i]))])
			
		if not label is None:
			hs,ls = GetLegendHandLab(ax)
			hs.append(ln[0])
			ls.append(label)
			ax.legend(hs,ls)
		
		ax.set_ylabel('$z_{SIII}$ (R$_J$)')
		ax.set_xlabel('$x_{SIII}$ (R$_J$)')


		mx = 1.1*mx	
		ax.set_xlim(-mx,mx)
		ax.set_ylim(-mx,mx)
		
		PlotJupiterXZ(ax)
		ax.set_aspect(1.0)

		return ax
	
	def PlotXY(self,ind='all',fig=None,maps=[1,1,0,0],label=None,color='black'):
		'''
		Plot field lines in the X-Y plane
		
		Inputs
		======
		ind : int|str
			Index of trace to plot. Can be scalar or an array. If set 
			ind='all' then all traces will be plotted.
		fig : None|pyplot|pyplot.Axes instance
			None - new figure will be created
			pyplot - new subplot will be created on existing figure
			pyplot.Axes - existing subplot will be used
		maps : list
			4-element array-like to determine the subplot position,
			ignored when fig=pyplot.Axes.
			maps = [xmaps,ymaps,xmap,ymap]
			xmaps - number of subplots in x-direction
			ymaps - number of subplots in y-direction
			xmap - x position of this subplot
			ymap - y position of this subplot
		label : None|str
			Add label to traces.
		color : str|array-like
			Colour to plot the field lines		
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
		
		x = self.x[ind]
		y = self.y[ind]

			
		mx = 1.5
		for i in range(0,ind.size):
			ln = ax.plot(y[i],x[i],color=color)
			mx = np.nanmax([mx,np.nanmax(np.abs(x[i])),np.nanmax(np.abs(y[i]))])
		if not label is None:
			hs,ls = GetLegendHandLab(ax)
			hs.append(ln[0])
			ls.append(label)
			ax.legend(hs,ls)
		yl = ax.get_xlim()
		ax.set_xlim(yl[::-1])
		
		ax.set_xlabel('$y_{SIII}$ (R$_J$)')
		ax.set_ylabel('$x_{SIII}$ (R$_J$)')

		mx = 1.1*mx	
		ax.set_xlim(mx,-mx)
		ax.set_ylim(-mx,mx)
		
		PlotJupiterXY(ax)
		ax.set_aspect(1.0)
		return ax
	
	def PlotRhoZ(self,ind='all',fig=None,maps=[1,1,0,0],label=None,color='black'):
		'''
		Plot field lines in the rho-Z plane

		
		Inputs
		======
		ind : int|str
			Index of trace to plot. Can be scalar or an array. If set 
			ind='all' then all traces will be plotted.
		fig : None|pyplot|pyplot.Axes instance
			None - new figure will be created
			pyplot - new subplot will be created on existing figure
			pyplot.Axes - existing subplot will be used
		maps : list
			4-element array-like to determine the subplot position,
			ignored when fig=pyplot.Axes.
			maps = [xmaps,ymaps,xmap,ymap]
			xmaps - number of subplots in x-direction
			ymaps - number of subplots in y-direction
			xmap - x position of this subplot
			ymap - y position of this subplot
		label : None|str
			Add label to traces.
		color : str|array-like
			Colour to plot the field lines		
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
		x = self.x[ind]
		y = self.y[ind]
		z = self.z[ind]

		
		r = np.array([np.sqrt(x[i]**2 + y[i]**2) for i in range(0,x.shape[0])],dtype='object')
		mx = 1.5
		for i in range(0,ind.size):
			ln = ax.plot(r[i],z[i],color=color)
			mx = np.nanmax([mx,np.nanmax(np.abs(r[i])),np.nanmax(np.abs(z[i]))])
		if not label is None:
			hs,ls = GetLegendHandLab(ax)
			hs.append(ln[0])
			ls.append(label)
			ax.legend(hs,ls)
		
		ax.set_ylabel('$z_{SIII}$ (R$_J$)')
		ax.set_xlabel(r'$\rho_{SIII}$ (R$_J$)')

		mx = 1.1*mx				
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
			for a in AI:
				ax.plot(self.s[t],self.halpha[t,a],label=r'Trace {:d} $\alpha=${:5.1f}'.format(t,self.alpha[a]))

		ax.legend()
		ax.set_xlabel(r'$s$ (R$_J$)')
		ax.set_ylabel(r'$h_{\alpha}$')

		return ax
		
		
	def _HourPos(self,lon,lat):
		'''
		Interpolate the positions each hour in the time range.
		
		'''
		if not self.Time:
			return None,None,None
		
		#get hours
		utch = np.unique(np.int32(self.utc)).astype('float64')
		use = np.where((utch >= self.utc[0]) & (utch <= self.utc[-1]))[0]
		utch = utch[use]
				
			
		#convert to Cartesian
		r = 90 - np.abs(lat)
		t = lon*np.pi/180.0
		x = r*np.cos(t)
		y = r*np.sin(t)
				
		#create interpolation objects
		fx = interp1d(self.utc,x)
		fy = interp1d(self.utc,y)
		
		#work out the position at those times
		xh = fx(utch)
		yh = fy(utch)
		
		rh = 90 - np.sqrt(xh**2 + yh**2)
		th = np.arctan2(yh,xh)
		
	
		return utch % 24,th,rh
		
	def StoreTime(self,**kwargs):
		'''
		Store time on the TraceField object (used for plotting)
		
		'''
		self.Time = False
		if 'Time' in kwargs:
			self._StoreTime(kwargs['Time'])
		elif 'Date' in kwargs and 'ut' in kwargs:
			self._StoreTime((kwargs['Date'],kwargs['ut']))

		
	def _StoreTime(self,Time):
		'''
		Store the time array in this object.
		
		'''
		if Time is None:
			#do nothing
			pass
		elif len(Time) == 2:
			#given Date and ut
			self.Date,self.ut = Time
			self.utc = TT.ContUT(self.Date,self.ut)
			self.Time = True
		else:
			#assume continuous time
			self.utc = Time
			self.Date,self.ut = TT.ContUTtoDate(self.utc)
			self.Time = True
		
	def PlotPigtail(self,Proj='normal',ShowLabels=True,Time=None,
					Date=None,ut=None,
					Hemisphere='both',colatlim=None,
					fig=None,maps=[1,1,0,0],**kwargs):
		'''
		Pigtail plot. I don't think it works.
		
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
		colatlim : None|float
			Limit of colatitude on the plot
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
			
		if not Time is None:
			self._StoreTime(Time)
		elif not Date is None and not ut is None:
			self._StoreTime((Date,ut))
		
		if fig is None:
			fig = plt
			fig.figure()
		if hasattr(fig,'Axes'):	
			ax = fig.subplot2grid((maps[1],maps[0]),(maps[3],maps[2]),projection='polar')
		else:
			ax = fig		
		ax.set_rlabel_position(0.0)
		rpo = ax.get_rlabel_position()
		ax.set_theta_zero_location("N")
		ax.set_rlim(rlim)
		if Hemisphere.lower() in ['both','north']:
			ax.plot(tn,rn,linewidth=kwargs.get('linewidth',2.0),color=kwargs.get('color','red'),label='North')	
			lono,lato = JupiterOvalNorth()
			ax.plot(lono*np.pi/180.0,lato,color='black',linestyle=':')
			if self.Time:
				uth,th,rh = self._HourPos(self.LonN,self.LatN)
				ax.scatter(th,rh,color=kwargs.get('color','red'),marker='o')
				for i in range(0,uth.size):
					ax.text(th[i],rh[i],'{:02d}'.format(np.int32(uth[i])),va='bottom',ha='left',color=kwargs.get('color','red'))
		if Hemisphere.lower() in ['both','south']:
			ax.plot(ts,rs,linewidth=kwargs.get('linewidth',2.0),color=kwargs.get('color','orange'),label='South')	
			lono,lato = JupiterOvalSouth()
			ax.plot(lono*np.pi/180.0,-lato,color='black',linestyle=':')
			if self.Time:
				uth,th,rh = self._HourPos(self.LonS,self.LatS)
				ax.scatter(th,rh,color=kwargs.get('color','orange'),marker='o')
				for i in range(0,uth.size):
					ax.text(th[i],rh[i],'{:02d}'.format(np.int32(uth[i])),va='bottom',ha='left',color=kwargs.get('color','orange'))
		ax.text((rpo+10.0)*np.pi/180.0,np.mean(rlim),'Latitude ($^\circ$)',rotation=rpo+90.0,ha='center',va='center')
		ax.set_xlabel('Longitude ($^\circ$)')
		ax.legend()
		
		return ax
		
	def TraceDict(self,RemoveNAN=True):
		'''
		Return a dictionary with all of the outputs of the field trace.
		
		Inputs
		======
		RemoveNAN : bool
			If True then arrays will be shortened by removing nans.
			
		Returns
		=======
		out : dict
			Contains the field traces coordinates, field components etc.
		
		'''
		#we could save a fair bit of space by removing NANs - this will
		#mean that simple 2D arrays will become arrays of objects
		if RemoveNAN:
			ptrs = ['x','y','z','Bx','By','Bz','s','R','Rnorm']
			out = {}
			keys = list(self.__dict__.keys())
			for k in keys:
				if k in ptrs:
					#2D
					tmp = np.zeros(self.n,dtype='object')
					for i in range(0,self.n):
						tmp[i] = self.__dict__[k][i,:self.nstep[i]]
					out[k] = tmp
				elif k == 'halpha':
					#3D
					tmp = np.zeros(self.halpha.shape[:2],dtype='object')
					for i in range(0,self.n):
						for j in range(0,self.nalpha):
							tmp[i,j] = self.halpha[i,j,:self.nstep[i]]
					out[k] = tmp
				else:
					out[k] = self.__dict__[k]
		else:
			out = self.__dict__
			
		#remove ctypes references
		out.pop('IntModelCode')
		out.pop('ExtModelCode')
			
		return out
	
	def Save(self,fname,RemoveNAN=True):
		'''
		Save the data in this object to file.
		
		Inputs
		======
		fname : str
			Path to the file where this trace will be save on disk.
		RemoveNAN : bool
			If True then arrays will be shortened by removing nans.
			
		'''
		out = self.TraceDict(RemoveNAN)
		
		print('Saving file: {:s}'.format(fname))
		
		pf.SaveObject(out,fname)


	def GetTrace(self,i):
		'''
		Return a trace.
		
		Inputs
		======
		i : int
			Index of the trace to be returned.

		
		Returns
		=======
		x : float
			x-coordinate (R_j)
		y : float
			y-coordinate (R_j)
		z : float
			z-coordinate (R_j)
		bx : float
			x-component of the magnetic field (nT)
		by : float
			y-component of the magnetic field (nT)
		bz : float
			z-component of the magnetic field (nT)
		r : float
			radial distance (R_E)
		rnorm : float
			Normalised radial distance (Rnorm = 1.0 at Rmax)
		s : float
			Distance along the field line trace (R_j)
		h : float
			H_alpha array.
		
		'''

		x = self.x[i][:self.nstep[i]]
		y = self.y[i][:self.nstep[i]]
		z = self.z[i][:self.nstep[i]]
		bx = self.Bx[i][:self.nstep[i]]
		by = self.By[i][:self.nstep[i]]
		bz = self.Bz[i][:self.nstep[i]]
		r = self.R[i][:self.nstep[i]]
		rnorm = self.Rnorm[i][:self.nstep[i]]
		s = self.s[i][:self.nstep[i]]
		h = self.halpha[i,:][:self.nstep[i]]
			
		return (x,y,z,bx,by,bz,r,rnorm,s,h)
