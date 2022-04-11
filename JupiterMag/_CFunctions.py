import numpy as np
import os
from . import Globals
from .ct import c_char_p,c_char_p_ptr
from .ct import c_bool,c_bool_ptr
from .ct import c_int,c_int_ptr
from .ct import c_float,c_float_ptr
from .ct import c_double,c_double_ptr,c_double_ptr_ptr
from ._CppLib import _GetLib

libjupitermag = _GetLib()

#field tracing routine
_CTraceField = libjupitermag.TraceField
_CTraceField.restype = None
_CTraceField.argtypes = [	c_int,				#n
							c_double_ptr,		#x0
							c_double_ptr,		#y0
							c_double_ptr,		#z0
							c_char_p,			#IntFunc
							c_int,				#Number of external functions
							c_char_p_ptr,		#ExtFunc
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
