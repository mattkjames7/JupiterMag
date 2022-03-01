import numpy as np
import ctypes as ct
import os
from .. import Globals
from .._CFunctions import libjupitermag

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
_CInternalField = libjupitermag.InternalField
_CInternalField.restype = None
_CInternalField.argtypes = [	c_int,			#number of elements
								c_double_ptr,	#x/r array
								c_double_ptr,	#y/t array
								c_double_ptr,	#z/p array
								c_double_ptr,	#Bx/Br output array
								c_double_ptr,	#Bx/Br output array
								c_double_ptr]	#Bx/Br output array
								
#internal model wrapper function
_CInternalFieldDeg = libjupitermag.InternalFieldDeg
_CInternalFieldDeg.restype = None
_CInternalFieldDeg.argtypes = [	c_int,			#number of elements
								c_double_ptr,	#x/r array
								c_double_ptr,	#y/t array
								c_double_ptr,	#z/p array
								c_int,			#MaxDeg
								c_double_ptr,	#Bx/Br output array
								c_double_ptr,	#Bx/Br output array
								c_double_ptr]	#Bx/Br output array


#set internal config
_CSetInternalCFG = libjupitermag.SetInternalCFG
_CSetInternalCFG.restype = None
_CSetInternalCFG.argtypes = [	c_char_p,
								c_bool,
								c_bool, 	
								c_int ]

#get internal config
_CGetInternalCFG = libjupitermag.GetInternalCFG
_CGetInternalCFG.restype = None
_CGetInternalCFG.argtypes = [	c_char_p,
								c_bool_ptr,
								c_bool_ptr,
								c_int_ptr 	]
