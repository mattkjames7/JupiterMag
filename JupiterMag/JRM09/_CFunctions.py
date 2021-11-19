import numpy as np
import ctypes as ct
import os
from .. import Globals
from .. _CFunctions import libinternal

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

#model field
_CJRM09FieldArray = libinternal.JRM09FieldArray
_CJRM09FieldArray.restype = None
_CJRM09FieldArray.argtypes = [	c_int, 			#number of elements
								c_double_ptr,	#p0
								c_double_ptr,	#p1
								c_double_ptr,	#p2
								c_double_ptr,	#B0
								c_double_ptr,	#B1
								c_double_ptr]	#B2

#set model config
_CSetJRM09Config = libinternal.SetJRM09Config
_CSetJRM09Config.restype = None
_CSetJRM09Config.argtypes = [	c_bool,			#Cart In
								c_bool]			#Cart Out
								
#get model config
_CGetJRM09Config = libinternal.GetJRM09Config
_CGetJRM09Config.restype = None
_CGetJRM09Config.argtypes = [	c_bool_ptr,		#Cart In
								c_bool_ptr ]	#Cart Out
