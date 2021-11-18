#ifndef __LIBCON2020_H__
#define __LIBCON2020_H__
#include <stdio.h>
#include <stdlib.h>
#include "con2020.h"
#include <string.h>

#endif
using namespace std;

/* we want to initialize the model objects with its parameters */
Con2020 con2020;

extern "C" {
	/* these wrappers can be used to get the magnetic field vectors */
	void Con2020Field(int n, double *p0, double *p1, double *p2,
					double *B0, double *B1, double *B2);
	
	void Con2020FieldVector(double p0, double p1, double p2,
			double *B0, double *B1, double *B2);


	void GetCon2020Params(double *mui, double *irho, double *r0, double *r1,
					double *d, double *xt, double *xp, char *eqtype,
					bool *Edwards, bool *ErrChk, bool *CartIn, bool *CartOut);
						
	
	void SetCon2020Params(double mui, double irho, double r0, double r1,
					double d, double xt, double xp, const char *eqtype,
					bool Edwards, bool ErrChk, bool CartIn, bool CartOut);

}
