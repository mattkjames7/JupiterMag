#include "libcon2020.h"

void Con2020Field(int n, double *p0, double *p1, double *p2,
			double *B0, double *B1, double *B2) {

	/* could create a separate function for default model */
	con2020.Field(n,p0,p1,p2,B0,B1,B2);
						
}

void Con2020FieldVector(double p0, double p1, double p2,
			double *B0, double *B1, double *B2) {

	/* could create a separate function for default model */
	con2020.Field(p0,p1,p2,B0,B1,B2);
						
}

void GetCon2020Params(double *mui, double *irho, double *r0, double *r1,
				double *d, double *xt, double *xp, char *eqtype,
				bool *Edwards, bool *ErrChk, bool *CartIn, bool *CartOut) {
	
	mui[0] = con2020.GetCurrentDensity();
	irho[0] = con2020.GetRadCurrentDensity();
	r0[0] = con2020.GetR0();
	r1[0] = con2020.GetR1();
	d[0] = con2020.GetCSHalfThickness();
	xt[0] = con2020.GetCSTilt();
	xp[0] = con2020.GetCSTiltAzimuth();
	Edwards[0] = con2020.GetEdwardsEqs();
	ErrChk[0] = con2020.GetErrCheck();
	CartIn[0] = con2020.GetCartIn();
	CartOut[0] = con2020.GetCartOut();
	con2020.GetEqType(eqtype);	
}
void SetCon2020Params(double mui, double irho, double r0, double r1,
				double d, double xt, double xp, const char *eqtype,
				bool Edwards, bool ErrChk, bool CartIn, bool CartOut) {
	
	con2020.SetCurrentDensity(mui);
	con2020.SetRadCurrentDensity(irho);
	con2020.SetR0(r0);
	con2020.SetR1(r1);
	con2020.SetCSHalfThickness(d);
	con2020.SetCSTilt(xt);
	con2020.SetCSTiltAzimuth(xp);
	con2020.SetEdwardsEqs(Edwards);
	con2020.SetErrCheck(ErrChk);
	con2020.SetCartIn(CartIn);
	con2020.SetCartOut(CartOut);
	con2020.SetEqType(eqtype);
}
