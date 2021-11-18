#include "libinternal.h"

void InternalField(int l, double *p0, double *p1, double *p2,
					double *B0, double *B1, double *B2,
					bool PolIn, bool PolOut, const char *Model) {

	/* select the appropriate model (default to VIP4) */
	if (strcmp(Model,"VIP4") == 0) {
		vip4.Field(l,p0,p1,p2,B0,B1,B2);
	} else if (strcmp(Model,"JRM09") == 0) {
		jrm09.Field(l,p0,p1,p2,B0,B1,B2);
	} else {
		vip4.Field(l,p0,p1,p2,B0,B1,B2);
	}		
						
}


void JRM09FieldArray(	int n, double *p0, double *p1, double *p2,
					double *B0, double *B1, double *B2) {
	jrm09.Field(n,p0,p1,p2,B0,B1,B2);
}

void JRM09Field(	double p0, double p1, double p2,
					double *B0, double *B1, double *B2) {
						
	jrm09.Field(p0,p1,p2,B0,B1,B2);
}

void SetJRM09Config(bool CartIn, bool CartOut) {
	
	jrm09.SetCartIn(CartIn);
	jrm09.SetCartOut(CartOut);
}

void GetJRM09Config(bool *CartIn, bool *CartOut) {
	
	CartIn[0] = jrm09.GetCartIn();
	CartOut[0] = jrm09.GetCartOut();
}

void VIP4FieldArray(	int n, double *p0, double *p1, double *p2,
					double *B0, double *B1, double *B2) {
	vip4.Field(n,p0,p1,p2,B0,B1,B2);
}

void VIP4Field(	double p0, double p1, double p2,
					double *B0, double *B1, double *B2) {
						
	vip4.Field(p0,p1,p2,B0,B1,B2);
}

void SetVIP4Config(bool CartIn, bool CartOut) {
	
	vip4.SetCartIn(CartIn);
	vip4.SetCartOut(CartOut);
}

void GetVIP4Config(bool *CartIn, bool *CartOut) {
	
	CartIn[0] = vip4.GetCartIn();
	CartOut[0] = vip4.GetCartOut();
}
