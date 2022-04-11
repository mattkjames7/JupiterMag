import numpy as np
import os
from .. import Globals
from .._CFunctions import libjupitermag

from ..ct import c_char_p,c_char_p_ptr
from ..ct import c_bool,c_bool_ptr
from ..ct import c_int,c_int_ptr
from ..ct import c_float,c_float_ptr
from ..ct import c_double,c_double_ptr,c_double_ptr_ptr

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
