import numpy as np


def _ptr2D(x):
	
	return (x.__array_interface__['data'][0] + np.arange(x.shape[0])*x.strides[0]).astype(np.uintp)
