import numpy as np
import os
from .. import Globals
import PyFileIO as pf

def _ReadTestPos():
	
	fname = Globals.ModulePath+"__data/testpos.dat"
	return pf.ReadASCIIData(fname)
