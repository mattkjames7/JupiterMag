#include "libinternal.h"

/* we want to initialize the model objects witht heir parameters */

InternalModel internalModel;

void InternalField(int n, double *p0, double *p1, double *p2,
					double *B0, double *B1, double *B2) {

	internalModel.Field(n,p0,p1,p2,B0,B1,B2);
				
}

void InternalFieldDeg(int n, double *p0, double *p1, double *p2,
					int MaxDeg, double *B0, double *B1, double *B2) {

	internalModel.Field(n,p0,p1,p2,MaxDeg,B0,B1,B2);
				
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

void GSFC13EVField(	double p0, double p1, double p2,
					double *B0, double *B1, double *B2) {
						
	gsfc13ev.Field(p0,p1,p2,B0,B1,B2);
}

void GSFC15EVField(	double p0, double p1, double p2,
					double *B0, double *B1, double *B2) {
						
	gsfc15ev.Field(p0,p1,p2,B0,B1,B2);
}

void GSFC15EVSField(	double p0, double p1, double p2,
					double *B0, double *B1, double *B2) {
						
	gsfc15evs.Field(p0,p1,p2,B0,B1,B2);
}

void ISAACField(	double p0, double p1, double p2,
					double *B0, double *B1, double *B2) {
						
	isaac.Field(p0,p1,p2,B0,B1,B2);
}

void JPL15EVField(	double p0, double p1, double p2,
					double *B0, double *B1, double *B2) {
						
	jpl15ev.Field(p0,p1,p2,B0,B1,B2);
}

void JPL15EVSField(	double p0, double p1, double p2,
					double *B0, double *B1, double *B2) {
						
	jpl15evs.Field(p0,p1,p2,B0,B1,B2);
}

void O4Field(	double p0, double p1, double p2,
					double *B0, double *B1, double *B2) {
						
	o4.Field(p0,p1,p2,B0,B1,B2);
}

void O6Field(	double p0, double p1, double p2,
					double *B0, double *B1, double *B2) {
						
	o6.Field(p0,p1,p2,B0,B1,B2);
}

void P11AField(	double p0, double p1, double p2,
					double *B0, double *B1, double *B2) {
						
	p11a.Field(p0,p1,p2,B0,B1,B2);
}

void SHAField(	double p0, double p1, double p2,
					double *B0, double *B1, double *B2) {
						
	sha.Field(p0,p1,p2,B0,B1,B2);
}

void U17EVField(	double p0, double p1, double p2,
					double *B0, double *B1, double *B2) {
						
	u17ev.Field(p0,p1,p2,B0,B1,B2);
}

void V117EVField(	double p0, double p1, double p2,
					double *B0, double *B1, double *B2) {
						
	v117ev.Field(p0,p1,p2,B0,B1,B2);
}

void VIPALField(	double p0, double p1, double p2,
					double *B0, double *B1, double *B2) {
						
	vipal.Field(p0,p1,p2,B0,B1,B2);
}

void VIT4Field(	double p0, double p1, double p2,
					double *B0, double *B1, double *B2) {
						
	vit4.Field(p0,p1,p2,B0,B1,B2);
}

