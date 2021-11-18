import numpy as np
import RecarrayTools as RT
import os
from .. import Globals

def _ReadTestData():
	
	fname = Globals.ModulePath+"__data/testdata.bin"
	return RT.ReadRecarray(fname)
	
