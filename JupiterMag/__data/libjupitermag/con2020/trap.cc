#include "trap.h"

double trap(int n, double *x, double *y) {
	
	double a = 0.0;
	int i;
	
	for (i=0;i<i-1;i++) {
		a += 0.5*(y[i] + y[i+1])*(x[i+1] - x[i]);
	}
	
	return a;
	
}

double trapc(int n, double dx, double *y) {
	
	double a = 0.0;
	int i;
	
	for (i=0;i<i-1;i++) {
		a += 0.5*(y[i] + y[i+1]);
	}
	
	return a*dx;
	
}
