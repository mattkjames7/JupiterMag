#ifndef __LIBINTERNAL_H__
#define __LIBINTERNAL_H__
#include <stdio.h>
#include <stdlib.h>
#include "internal.h"
#include <string.h>

#endif
using namespace std;

/* we want to initialize the model objects witht heir parameters */
Internal vip4("VIP4");
Internal jrm09("JRM09");

extern "C" {
	/* these wrappers can be used to get the magnetic field vectors */
	void InternalField(int l, double *p0, double *p1, double *p2,
					double *B0, double *B1, double *B2,
					bool PolIn, bool PolOut, const char *Model);


	void JRM09FieldArray(	int n, double *p0, double *p1, double *p2,
						double *B0, double *B1, double *B2);
	void JRM09Field(	double p0, double p1, double p2,
						double *B0, double *B1, double *B2);
	void SetJRM09Config(bool CartIn, bool CartOut);
	void GetJRM09Config(bool *CartIn, bool *CartOut);
	void VIP4FieldArray(	int n, double *p0, double *p1, double *p2,
						double *B0, double *B1, double *B2);
	void VIP4Field(	double p0, double p1, double p2,
						double *B0, double *B1, double *B2);
	void SetVIP4Config(bool CartIn, bool CartOut);
	void GetVIP4Config(bool *CartIn, bool *CartOut);
}
