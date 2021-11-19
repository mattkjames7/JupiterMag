import numpy as np
from ._CFunctions import _CTraceField
import ctypes
from ._ptr2D import _ptr2D
from . import JRM09
from . import VIP4 
from . import Con2020
		
class TraceField(object):
	'''
	Object which stores the result of a magnetic field trace or a series 
	of traces performed using a combination of internal and external field
	models
	
	
	'''
	
	def __init__(self,x0,y0,z0,IntModel='JRM09',ExtModel='Con2020', 
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
		self.x0 = np.array(x0).astype("float64")
		self.y0 = np.array(y0).astype("float64")
		self.z0 = np.array(z0).astype("float64")
		self.n = np.int32(self.x0.size)

		self.IntModel = IntModel
		self.IntModelCode = ctypes.c_char_p(IntModel.encode('utf-8'))
		self.ExtModel = ExtModel
		self.ExtModelCode = ctypes.c_char_p(ExtModel.encode('utf-8'))

		#make sure models are in Cartesian
		Models = [IntModel,ExtModel]
		if "JRM09" in Models:
			JRM09.Config(CartesianIn=True,CartesianOut=True)
		elif "VIP4" in Models:
			VIP4.Config(CartesianIn=True,CartesianOut=True)
		elif "Con2020" in Models:
			Con2020.Config(CartesianIn=True,CartesianOut=True)
			

		#kwargs
		defargs = {	'MaxLen'	:		1000,
					'MaxStep'	:		1.0,
					'InitStep'	:		0.5,
					'MinStep'	:		0.001,
					'ErrMax'	:		0.0001,
					'Delta'		:		0.05,
					'Verbose'	:		False,
					'TraceDir'	:		'both',
					'alpha' 	:		[0.0,90.0]}
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

		
	
		
	
