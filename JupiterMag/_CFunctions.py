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
	libinternal = ct.CDLL(os.path.dirname(__file__)+"/__data/libjupitermag/internal/libinternal.so")

# check that the con2020 field library exists
try:
	libcon2020 = ct.CDLL(os.path.dirname(__file__)+"/__data/libjupitermag/con2020/libcon2020.so")
except:
	print('importing libcon2020.so failed, attempting to recompile')
	path = os.path.dirname(__file__)
	if '/usr/local/' in path:
		sudo = 'sudo '
	else:
		sudo = ''

	CWD = os.getcwd()
	os.chdir(os.path.dirname(__file__)+"/__data/libjupitermag/con2020/")
	os.system(sudo+'make clean')
	os.system(sudo+'make')
	os.chdir(CWD)	
	libcon2020 = ct.CDLL(os.path.dirname(__file__)+"/__data/libjupitermag/con2020/libcon2020.so")

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


#con2020 wrapper function
_CCon2020Field = libcon2020.Con2020Field
_CCon2020Field.restype = None
_CCon2020Field.argtypes = [		c_int,			#number of input elements
								c_double_ptr,	#x/r array
								c_double_ptr,	#y/t array
								c_double_ptr,	#z/p array
								c_double_ptr,	#Bx/Br output array
								c_double_ptr,	#Bx/Br output array
								c_double_ptr]	#Bx/Br output array]
								
_CCon2020Vector = libcon2020.Con2020Vector
_CCon2020Vector.restype = None
_CCon2020Vector.argtypes = [	c_double,	#x/r scalar
								c_double,	#y/t scalar
								c_double,	#z/p scalar
								c_double_ptr,	#Bx/Br output scalar
								c_double_ptr,	#Bx/Br output scalar
								c_double_ptr]	#Bx/Br output scalar

#con2020 get parameters
_CGetCon2020Params = libcon2020.GetCon2020Params
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
								c_bool_ptr ]		#CartOut
#con2020 set parameters
_CSetCon2020Params = libcon2020.SetCon2020Params
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
								c_bool ]		#CartOut
