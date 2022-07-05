import time
from .. import Internal
import numpy as np



def Timer():
	
	ns = [1,10,100,1000,10000,100000]
	l = len(ns)
	
	mu_vec_pol = np.zeros(l,dtype='float64')
	sd_vec_pol = np.zeros(l,dtype='float64')
	mu_sgl_pol = np.zeros(l,dtype='float64')
	sd_sgl_pol = np.zeros(l,dtype='float64')
	mu_vec_car = np.zeros(l,dtype='float64')
	sd_vec_car = np.zeros(l,dtype='float64')
	mu_sgl_car = np.zeros(l,dtype='float64')
	sd_sgl_car = np.zeros(l,dtype='float64')
	
	for i in range(0,l):
		print('n = ',ns[i])
		mu_vec_pol[i],mu_sgl_pol[i],mu_vec_car[i],mu_sgl_car[i],sd_vec_pol[i],sd_sgl_pol[i],sd_vec_car[i],sd_sgl_car[i] = TimeN(ns[i])
		
	print('Mean Vectorized Polar: ',' '.join(mu_vec_pol.astype('U')))
	print('StdDev Vectorized Polar: ',' '.join(sd_vec_pol.astype('U')))
	print('Mean Single Polar: ',' '.join(mu_sgl_pol.astype('U')))
	print('StdDev Single Polar: ',' '.join(sd_sgl_pol.astype('U')))
		
	print('Mean Vectorized Cartesian: ',' '.join(mu_vec_car.astype('U')))
	print('StdDev Vectorized Cartesian: ',' '.join(sd_vec_car.astype('U')))
	print('Mean Single Cartesian: ',' '.join(mu_sgl_car.astype('U')))
	print('StdDev Single Cartesian: ',' '.join(sd_sgl_car.astype('U')))
	


def TimeN(n):
	'''
	time how long a model takes to calcualte vectors
	
	'''
	ntest = 5
	
	#get some randomized vectors
	x = np.random.rand(n)*100 - 50.0
	y = np.random.rand(n)*100 - 50.0
	z = np.random.rand(n)*40 - 20.0
	
	r = np.sqrt(x**2 + y**2 + z**2)
	t = np.arccos(z/r)
	p = np.arctan2(y,x)
	
	#storage
	dt_vec_pol = np.zeros(ntest,dtype='float64')
	dt_sgl_pol = np.zeros(ntest,dtype='float64')
	dt_vec_car = np.zeros(ntest,dtype='float64')
	dt_sgl_car = np.zeros(ntest,dtype='float64')
	
	#set model
	Model = 'jrm33'
	Degree =  18
	Internal.Config(Model=Model,Degree=Degree)
	print('Testing model {:s}, degree {:d}'.format(Model,Degree))
	
	#test vectorized polar
	Internal.Config(CartesianIn=False,CartesianOut=False)
	for i in range(0,ntest):
		t0 = time.time()
		B0,B1,B2 = Internal.Field(r,t,p)
		t1 = time.time()
		dt_vec_pol[i] = t1 - t0
	mu_vec_pol = np.mean(dt_vec_pol)
	sd_vec_pol = np.std(dt_vec_pol)
	print("{:12.10f} +/- {:12.10f} s - polar, vectorized ({:d} vectors, {:d} tests)".format(mu_vec_pol,sd_vec_pol,n,ntest))
	
	#test scalar polar
	Internal.Config(CartesianIn=False,CartesianOut=False)
	for i in range(0,ntest):
		t0 = time.time()
		for j in range(0,n):
			B0,B1,B2 = Internal.Field(r[j],t[j],p[j])
		t1 = time.time()
		dt_sgl_pol[i] = t1 - t0
	mu_sgl_pol = np.mean(dt_sgl_pol)
	sd_sgl_pol = np.std(dt_sgl_pol)
	print("{:12.10f} +/- {:12.10f} s - polar, single vectors ({:d} vectors, {:d} tests)".format(mu_sgl_pol,sd_sgl_pol,n,ntest))

	
	#test vectorized cart
	Internal.Config(CartesianIn=True,CartesianOut=True)
	for i in range(0,ntest):
		t0 = time.time()
		B0,B1,B2 = Internal.Field(x,y,z)
		t1 = time.time()
		dt_vec_car[i] = t1 - t0
	mu_vec_car = np.mean(dt_vec_car)
	sd_vec_car = np.std(dt_vec_car)
	print("{:12.10f} +/- {:12.10f} s - Cartesian, vectorized ({:d} vectors, {:d} tests)".format(mu_vec_car,sd_vec_car,n,ntest))
	
	#test scalar cart
	Internal.Config(CartesianIn=True,CartesianOut=True)
	for i in range(0,ntest):
		t0 = time.time()
		for j in range(0,n):
			B0,B1,B2 = Internal.Field(x[j],y[j],z[j])
		t1 = time.time()
		dt_sgl_car[i] = t1 - t0
	mu_sgl_car = np.mean(dt_sgl_car)
	sd_sgl_car = np.std(dt_sgl_car)
	print("{:12.10f} +/- {:12.10f} s - Cartesian, single vectors ({:d} vectors, {:d} tests)".format(mu_sgl_car,sd_sgl_car,n,ntest))
	
			
		
	return mu_vec_pol,mu_sgl_pol,mu_vec_car,mu_sgl_car,sd_vec_pol,sd_sgl_pol,sd_vec_car,sd_sgl_car
