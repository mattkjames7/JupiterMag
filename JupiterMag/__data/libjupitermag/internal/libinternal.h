#ifndef __LIBINTERNAL_H__
#define __LIBINTERNAL_H__
#include <stdio.h>
#include <stdlib.h>
#include "internal.h"
#include <string.h>
#include "internalmodel.h"

using namespace std;

/* we want to initialize the model objects witht heir parameters */
extern Internal vip4;
extern Internal jrm09;
extern InternalModel internalModel;

extern "C" {
	/* these wrappers can be used to get the magnetic field vectors */

	void InternalField(int n, double *p0, double *p1, double *p2,
						double *B0, double *B1, double *B2);
	void SetInternalCFG(char *Model, bool CartIn, bool CartOut);

	void GetInternalCFG(char *Model, bool *CartIn, bool *CartOut);

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
#endif
