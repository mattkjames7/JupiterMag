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
_CVIP4FieldArray = libinternal.VIP4FieldArray
_CVIP4FieldArray.restype = None
_CVIP4FieldArray.argtypes = [	c_int, 			#number of elements
								c_double_ptr,	#p0
								c_double_ptr,	#p1
								c_double_ptr,	#p2
								c_double_ptr,	#B0
								c_double_ptr,	#B1
								c_double_ptr]	#B2

#set model config
_CSetVIP4Config = libinternal.SetVIP4Config
_CSetVIP4Config.restype = None
_CSetVIP4Config.argtypes = [	c_bool,			#Cart In
								c_bool]			#Cart Out
								
#get model config
_CGetVIP4Config = libinternal.GetVIP4Config
_CGetVIP4Config.restype = None
_CGetVIP4Config.argtypes = [	c_bool_ptr,		#Cart In
								c_bool_ptr ]	#Cart Out
