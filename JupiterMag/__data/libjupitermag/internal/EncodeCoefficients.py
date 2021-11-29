import numpy as np
import os

def EncodeCoefficients(fnamein,fnameout):
	'''
	This will encode the ASCII files of internal magnetic field model
	coefficients as pure binary.
	
	Inputs
	======
	fnamein : str
		Name of the ASCII file containing the coefficients.
	fnameout : str
		Name of the output binary file.
	
	'''
	
	#open the ASCII file
	f = open(fnamein,'r')
	lines = f.readlines()
	f.close()
	lines = np.array(lines)
	
	#get the number of lines in the file
	nl = lines.size
	
	#create the arrays for the coefficients
	gh = np.zeros(nl,dtype='int8')
	n = np.zeros(nl,dtype='int32')
	m = np.zeros(nl,dtype='int32')
	coeff = np.zeros(nl,dtype='float64')
	
	#fill them
	for i in range(0,nl):
		s = lines[i].split()

		if s[0] == 'h':
			gh[i] = 1
		else:
			gh[i] = 0
			
		n[i] = np.int32(s[1])
		m[i] = np.int32(s[2])
		coeff[i] = np.float64(s[3])
	
	#open the output file
	print('Saving {:s}'.format(fnameout))
	f = open(fnameout,'wb')
	np.int32(nl).tofile(f)
	gh.tofile(f)
	n.tofile(f)
	m.tofile(f)
	coeff.tofile(f)
	f.close()
	

if __name__ == "__main__":
	coeffs = ['vipalcoeffs','isaaccoeffs','jrm09coeffs','vip4coeffs']
	for c in coeffs:
		print('Converting {:s}'.format(c))
		EncodeCoefficients(c+'.dat',c+'.bin')
