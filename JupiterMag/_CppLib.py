import numpy as np
import os
import subprocess
import ctypes
import platform

def _LibPath():
	'''
	Return a path to the C++ library
	
	Returns
	=======
	path : str
		path to the library file.
	
	'''
	return os.path.dirname(__file__)+"/__data/libjupitermag/lib/"

def _LibName(WithPath=False):
	'''
	Return the name of the library.
	
	Inputs
	======
	WithPath : bool
		If True then the fuull path to the library will be included.
		
	Returns
	=======
	libpath : str
		Library name.
	
	'''
	if WithPath:
		path = _LibPath()
	else:
		path = ''

	osname = platform.uname().system
	libexts = {	'Linux':'so',
				'Windows':'dll',
				'Darwin':'dylib'}	
	
	ext = libexts[osname]

	if ext is None:
		raise Exception("The Operating System ({:s}) is not supported".format(osname))
	
	return path + 'libjupitermag.' + ext


def _LibExists():
	'''
	Check if the library file exists.
	
	Returns
	=======
	exists : bool
		True if the file exists
	'''
	return os.path.isfile(_LibName(True))
	
def _GetLib():
	'''	
	Return an instance of the C++ library
	
	Returns
	=======
	lib : ctypes.CDLL
		C++ library containing the field model code
	'''
	fname = _LibName(True)
	
	try:
		print('Importing Library')
		lib = ctypes.CDLL(fname)
		print('done')
	except:
		print("Importing C++ library failed. Please reinstall...")
		raise SystemExit
		
	return lib

