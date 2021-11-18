#include "model.h"

void ModelField(double p0, double p1, double p2, 
				const char *internal, const char *external, 
				bool CartIn, bool CartOut,
				double *B0, double *B1, double *B2) {

	/* get the internal field model */
	double Bi0, Bi1, Bi2;
	if (strcmp(internal,"JRM09") == 0) {
		jrm09.SetCartIn(CartIn);
		jrm09.SetCartOut(CartOut);
		jrm09.Field(p0,p1,p2,&Bi0,&Bi1,&Bi2);
	} else if (strcmp(internal,"VIP4") == 0) {
		vip4.SetCartIn(CartIn);
		vip4.SetCartOut(CartOut);
		vip4.Field(p0,p1,p2,&Bi0,&Bi1,&Bi2);
	} else {
		Bi0 = 0.0;
		Bi1 = 0.0;
		Bi2 = 0.0;
	}
	
	/* and the external field */
	double Be0, Be1, Be2;
	if (strcmp(external,"Con2020") == 0) {
		con2020.SetCartIn(CartIn);
		con2020.SetCartOut(CartOut);
		con2020.Field(p0,p1,p2,&Be0,&Be1,&Be2);
	} else {
		Be0 = 0.0;
		Be1 = 0.0;
		Be2 = 0.0;
	}
	
	B0[0] = Bi0 + Be0;
	B1[0] = Bi1 + Be1;
	B2[0] = Bi2 + Be2;

}

void ModelFieldArray(	int n, double *p0, double *p1, double *p2, 
						const char *internal, const char *external, 
						bool CartIn, bool CartOut,
						double *B0, double *B1, double *B2) {

	/* get the internal field model */
	int i;
	double *Bi0 = new double[n];
	double *Bi1 = new double[n];
	double *Bi2 = new double[n];
	if (strcmp(internal,"JRM09") == 0) {
		jrm09.SetCartIn(CartIn);
		jrm09.SetCartOut(CartOut);
		jrm09.Field(n,p0,p1,p2,Bi0,Bi1,Bi2);
	} else if (strcmp(internal,"VIP4") == 0) {
		vip4.SetCartIn(CartIn);
		vip4.SetCartOut(CartOut);
		vip4.Field(n,p0,p1,p2,Bi0,Bi1,Bi2);
	} else {
		for (i=0;i<n;i++) {
			Bi0[i] = 0.0;
			Bi1[i] = 0.0;
			Bi2[i] = 0.0;
		}
	}
	
	/* and the external field */
	double *Be0 = new double[n];
	double *Be1 = new double[n];
	double *Be2 = new double[n];
	if (strcmp(external,"Con2020") == 0) {
		con2020.SetCartIn(CartIn);
		con2020.SetCartOut(CartOut);
		con2020.Field(n,p0,p1,p2,Be0,Be1,Be2);
	} else {
		for (i=0;i<n;i++) {
			Be0[i] = 0.0;
			Be1[i] = 0.0;
			Be2[i] = 0.0;
		}
	}
	
	for (i=0;i<n;i++) {
		B0[i] = Bi0[i] + Be0[i];
		B1[i] = Bi1[i] + Be1[i];
		B2[i] = Bi2[i] + Be2[i];
	}
	
	delete[] Bi0;
	delete[] Bi1;
	delete[] Bi2;
	delete[] Be0;
	delete[] Be1;
	delete[] Be2;
}
