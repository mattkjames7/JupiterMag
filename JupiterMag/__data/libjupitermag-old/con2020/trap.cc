#include "trap.h"

double trap(int n, double *x, double *y) {
	
	double a = 0.0;
	double tmp;
	int i;
	
	for (i=0;i<n-1;i++) {
		tmp = 0.5*(y[i] + y[i+1])*(x[i+1] - x[i]);
		if (!isnan(tmp)) {
			a += tmp;
		}
	}
	
	return a;
	
}

double trapc(int n, double dx, double *y) {
	
	double a = 0.0;
	double tmp;
	int i;
	
	for (i=0;i<n-1;i++) {
		
		tmp = 0.5*(y[i] + y[i+1]);
		if (!isnan(tmp)) {
			a += tmp;
		}
	}
	return a*dx;
	
}
