#include "libinternal.h"

/* we want to initialize the model objects witht heir parameters */
Internal vip4("VIP4");
Internal jrm09("JRM09");
InternalModel internalModel;

void InternalField(int n, double *p0, double *p1, double *p2,
					double *B0, double *B1, double *B2) {

	internalModel.Field(n,p0,p1,p2,B0,B1,B2);
				
}

void SetInternalCFG(char *Model, bool CartIn, bool CartOut) {
	internalModel.SetCartIn(CartIn);
	internalModel.SetCartOut(CartOut);
	internalModel.SetModel(Model);
}

void GetInternalCFG(char *Model, bool *CartIn, bool *CartOut) {
	CartIn[0] = internalModel.GetCartIn();
	CartOut[0] = internalModel.GetCartOut();
	internalModel.GetModel(Model);
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
