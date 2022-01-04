#ifndef __LIBINTERNAL_H__
#define __LIBINTERNAL_H__
#include <stdio.h>
#include <stdlib.h>
#include "internal.h"
#include <string.h>
#include "internalmodel.h"

using namespace std;

/* we want to initialize the model objects witht heir parameters */
extern InternalModel internalModel;

extern "C" {
	/* these wrappers can be used to get the magnetic field vectors */

	void InternalField(int n, double *p0, double *p1, double *p2,
						double *B0, double *B1, double *B2);
	void InternalFieldDeg(int n, double *p0, double *p1, double *p2,
						int MaxDeg, double *B0, double *B1, double *B2);
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
	

	void GSFC13EVField(	double p0, double p1, double p2,
						double *B0, double *B1, double *B2);

	void GSFC15EVField(	double p0, double p1, double p2,
						double *B0, double *B1, double *B2);

	void GSFC15EVSField(	double p0, double p1, double p2,
						double *B0, double *B1, double *B2);

	void ISAACField(	double p0, double p1, double p2,
						double *B0, double *B1, double *B2);

	void JPL15EVField(	double p0, double p1, double p2,
						double *B0, double *B1, double *B2);

	void JPL15EVSField(	double p0, double p1, double p2,
						double *B0, double *B1, double *B2);

	void O4Field(	double p0, double p1, double p2,
						double *B0, double *B1, double *B2);

	void O6Field(	double p0, double p1, double p2,
						double *B0, double *B1, double *B2);

	void P11AField(	double p0, double p1, double p2,
						double *B0, double *B1, double *B2);

	void SHAField(	double p0, double p1, double p2,
						double *B0, double *B1, double *B2);

	void U17EVField(	double p0, double p1, double p2,
						double *B0, double *B1, double *B2);

	void V117EVField(	double p0, double p1, double p2,
						double *B0, double *B1, double *B2);

	void VIPALField(	double p0, double p1, double p2,
						double *B0, double *B1, double *B2);

	void VIT4Field(	double p0, double p1, double p2,
						double *B0, double *B1, double *B2);

	
}
#endif
