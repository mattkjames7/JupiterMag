import numpy as np
import ctypes as ct
import os
from . import Globals

# # check that the internal field library exists
# try:
	# libinternal = ct.CDLL(Globals.ModulePath+"__data/libjupitermag/internal/libinternal.so")
# except:
	# print('importing libinternal.so failed, attempting to recompile')
	# path = Globals.ModulePath
	# if '/usr/local/' in path:
		# sudo = 'sudo '
	# else:
		# sudo = ''

	# CWD = os.getcwd()
	# os.chdir(Globals.ModulePath+"__data/libjupitermag/internal/")
	# os.system(sudo+'make clean')
	# os.system(sudo+'make')
	# os.chdir(CWD)	
	# libinternal = ct.CDLL(Globals.ModulePath+"__data/libjupitermag/internal/libinternal.so")

# check that the jupiter mag field library exists
try:
	libjupitermag = ct.CDLL(Globals.ModulePath+"__data/libjupitermag/libjupitermag.so")
except:
	print('importing libjupitermag.so failed, attempting to recompile')
	path = Globals.ModulePath
	if '/usr/local/' in path:
		sudo = 'sudo '
	else:
		sudo = ''

	CWD = os.getcwd()
	os.chdir(Globals.ModulePath+"__data/libjupitermag/")
	os.system(sudo+'make clean')
	os.system(sudo+'make')
	os.chdir(CWD)	
	libjupitermag = ct.CDLL(Globals.ModulePath+"__data/libjupitermag/libjupitermag.so")

#define some dtypes
c_char_p = ct.c_char_p
c_bool = ct.c_bool
c_int = ct.c_int
c_float = ct.c_float
c_double = ct.c_double
c_float_ptr = np.ctypeslib.ndpointer(ct.c_float,flags="C_CONTIGUOUS")
c_int_ptr = np.ctypeslib.ndpointer(ct.c_int,flags="C_CONTIGUOUS")
c_bool_ptr = np.ctypeslib.ndpointer(ct.c_bool,flags="C_CONTIGUOUS")
#this one is a hack found at: https://stackoverflow.com/a/32138619/15482422
#it allows us to send None instead of an array which is treated as NULL
c_double_ptr_base = np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS")
def _from_param(cls, obj):
		if obj is None:
			return obj
		return c_double_ptr_base.from_param(obj)
c_double_ptr = type('c_double_ptr',(c_double_ptr_base,),{'from_param':classmethod(_from_param)})

c_double_ptr_ptr = np.ctypeslib.ndpointer(np.uintp,ndim=1,flags="C_CONTIGUOUS")

#internal model wrapper function
_CInternalField = libjupitermag.InternalField
_CInternalField.restype = None
_CInternalField.argtypes = [	c_int,			#number of elements
								c_double_ptr,	#x/r array
								c_double_ptr,	#y/t array
								c_double_ptr,	#z/p array
								c_double_ptr,	#Bx/Br output array
								c_double_ptr,	#Bx/Br output array
								c_double_ptr,	#Bx/Br output array
								c_char_p ]		#Model string


#field tracing routine
_CTraceField = libjupitermag.TraceField
_CTraceField.restype = None
_CTraceField.argtypes = [	c_int,				#n
							c_double_ptr,		#x0
							c_double_ptr,		#y0
							c_double_ptr,		#z0
							c_char_p,			#IntFunc
							c_char_p,			#ExtFunc
							c_int,				#MaxLen
							c_double,			#MaxStep
							c_double,			#InitStep
							c_double,			#MinStep
							c_double,			#ErrMax
							c_double,			#Delta
							c_bool,				#Verbose
							c_int,				#TraceDir
							c_int_ptr,			#nstep
							c_double_ptr_ptr,	#x
							c_double_ptr_ptr,	#y
							c_double_ptr_ptr,	#z
							c_double_ptr_ptr,	#Bx
							c_double_ptr_ptr,	#By
							c_double_ptr_ptr,	#Bz
							c_double_ptr_ptr,	#R
							c_double_ptr_ptr,	#S
							c_double_ptr_ptr,	#Rnorm
							c_double_ptr_ptr,	#FP
							c_int,				#nalpha
							c_double_ptr,		#alpha
							c_double_ptr]		#halpha
