import numpy as np
import os
from .. import Globals
from .._CFunctions import libjupitermag

from ..ct import c_char_p,c_char_p_ptr
from ..ct import c_bool,c_bool_ptr
from ..ct import c_int,c_int_ptr
from ..ct import c_float,c_float_ptr
from ..ct import c_double,c_double_ptr,c_double_ptr_ptr


#con2020 wrapper function
_CCon2020FieldArray = libjupitermag.Con2020FieldArray
_CCon2020FieldArray.restype = None
_CCon2020FieldArray.argtypes = [c_int,			#number of input elements
								c_double_ptr,	#x/r array
								c_double_ptr,	#y/t array
								c_double_ptr,	#z/p array
								c_double_ptr,	#Bx/Br output array
								c_double_ptr,	#Bx/Br output array
								c_double_ptr]	#Bx/Br output array]
								
_CCon2020Field = libjupitermag.Con2020Field
_CCon2020Field.restype = None
_CCon2020Field.argtypes = [		c_double,	#x/r scalar
								c_double,	#y/t scalar
								c_double,	#z/p scalar
								c_double_ptr,	#Bx/Br output scalar
								c_double_ptr,	#Bx/Br output scalar
								c_double_ptr]	#Bx/Br output scalar

#con2020 get parameters
_CGetCon2020Params = libjupitermag.GetCon2020Params
_CGetCon2020Params.restype = None
_CGetCon2020Params.argtypes = [	c_double_ptr,		#mui
								c_double_ptr,		#irho
								c_double_ptr,		#r0
								c_double_ptr,		#r1
								c_double_ptr,		#d
								c_double_ptr,		#xt
								c_double_ptr,		#xp
								c_char_p,		#eqtype
								c_bool_ptr,			#Edwards
								c_bool_ptr,			#ErrChk
								c_bool_ptr,			#CartIn
								c_bool_ptr,			#CartOut
								c_bool_ptr, 		#Smooth
								c_double_ptr,		#DeltaRho
								c_double_ptr,		#Deltaz,
								c_double_ptr,		#g
								c_char_p,			#azfunc
								c_double_ptr,		#wO_open
								c_double_ptr,		#wO_om
								c_double_ptr,		#thetamm
								c_double_ptr,		#dthetamm
								c_double_ptr,		#thetaoc
								c_double_ptr]		#dthetaoc
#con2020 set parameters
_CSetCon2020Params = libjupitermag.SetCon2020Params
_CSetCon2020Params.restype = None
_CSetCon2020Params.argtypes = [	c_double,		#mui
								c_double,		#irho
								c_double,		#r0
								c_double,		#r1
								c_double,		#d
								c_double,		#xt
								c_double,		#xp
								c_char_p,		#eqtype
								c_bool,			#Edwards
								c_bool,			#ErrChk
								c_bool,			#CartIn
								c_bool,			#CartOut
								c_bool, 		#Smooth
								c_double,	#DeltaRho
								c_double,	#DeltaZ
								c_double,		#g
								c_char_p,			#azfunc
								c_double,		#wO_open
								c_double,		#wO_om
								c_double,		#thetamm
								c_double,		#dthetamm
								c_double,		#thetaoc
								c_double]		#dthetaoc