import numpy as np
import ctypes as ct
import os

# check that the internal field library exists
try:
	libinternal = ct.CDLL(os.path.dirname(__file__)+"/__data/libjupitermag/internal/libinternal.so")
except:
	print('importing libinternal.so failed, attempting to recompile')
	path = os.path.dirname(__file__)
	if '/usr/local/' in path:
		sudo = 'sudo '
	else:
		sudo = ''

	CWD = os.getcwd()
	os.chdir(os.path.dirname(__file__)+"/__data/libjupitermag/internal/")
	os.system(sudo+'make clean')
	os.system(sudo+'make')
	os.chdir(CWD)	
	liblsmodel = ct.CDLL(os.path.dirname(__file__)+"/__data/libjupitermag/internal/libinternal.so")

#define some dtypes
c_char_p = ct.c_char_p
c_bool = ct.c_bool
c_int = ct.c_int
c_float = ct.c_float
c_double = ct.c_double
c_float_ptr = np.ctypeslib.ndpointer(ct.c_float,flags="C_CONTIGUOUS")
c_double_ptr = np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS")
c_int_ptr = np.ctypeslib.ndpointer(ct.c_int,flags="C_CONTIGUOUS")
c_bool_ptr = np.ctypeslib.ndpointer(ct.c_bool,flags="C_CONTIGUOUS")


#internal model wrapper function
_CInternalField = libinternal.InternalField
_CInternalField.restype = None
_CInternalField.argtypes = [	c_int,			#number of elements
								c_double_ptr,	#x/r array
								c_double_ptr,	#y/t array
								c_double_ptr,	#z/p array
								c_double_ptr,	#Bx/Br output array
								c_double_ptr,	#Bx/Br output array
								c_double_ptr,	#Bx/Br output array
								c_bool,			#Polar coordinate input
								c_bool,			#Polar coordinate output
								c_char_p ]		#Model string



