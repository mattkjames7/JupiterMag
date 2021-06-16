#include "cn2020.h"

Con2020::Con2020() {
	/* set all model parameters to their default values */
	mui_ = 139.6;
	irho_ = 16.7;
	r0_ = 7.8;
	r1_ = 51.4;
	d_ = 3.6;
	xt_ = 9.3;
	xp_ = -24.2;
	eqtype_ = "hybrid";
	Edwards_ = true;
	
	/* some other values which will only need calculating once */
	dipshift_ = xt_*deg2rad;
	diptilt_ = xp_*deg2rad;
	
	
}

void Con2020::_SysIII2Mag(int n, double *x0, double *y0, double *z0,
						double *x1, double *y1, double *z1) {
	
	
	/* some temporary variables which get used more than once */
	double sindt, cosdt, cosds, sinds, xt;
	sinds = sin(dipshift_);
	cosds = cos(dipshift_);
	sindt = sin(diptilt_);
	cosdt = cos(diptilt_);
	
	
	int i;
	for (i=0;i<n;i++) {
		/*intermediate value for x */
		xt = x0[i]*cosds + y0[i]*sinds[i];
		
		/*newly rotated coords */
		x1[i] = xt*cosdt + z[i]sindt;
		y1[i] = y[i]*cosds - x[i]*sinds;
		z1[i] = z[i]*cosdt - xt*sindt;
	}
	
					
}

void Con2020::_PolSysIII2Mag(int n, double *r, double *t, double *p,
						double *x1, double *y1, double *z1) {
	
	
	/* some temporary variables which get used more than once */
	double sindt, cosdt, cosp, sinp, rsint, rcost;

	sindt = sin(diptilt_);
	cosdt = cos(diptilt_);
	
	
	int i;
	for (i=0;i<n;i++) {
		/* temporary variables */
		rsint = r[i]*sin(t[i]);
		rcost = r[i]*cos(t[i]);
		sinp = sin(p[i] - dipshift_);
		cosp = cos(p[i] - dipshift_);		
		
		/*newly rotated coords */
		x1[i] = rsint*cosp*cosdt + rcost*sindt;
		y1[i] = rsint*sinp;
		z1[i] = rcost*rcosdt - rsint*cosp*sindt;
	}
	
					
}


void Con2020::Field(int n, double *p0, double *p1, double *p2, 
			double *B0, double *B1, double *B2,
			bool PolIn, bool PolOut) {

	/* allocate arrays to store model coordinates */
	double *x = new double[n];
	double *y = new double[n];
	double *z = new double[n];


	/* convert input coordinates to the coordinate system used by the model */
	if (PolIn) {
		/* this will convert spherical polar to the Cartesian model coordinates */
		_PolSysIII2Mag(n,p0,p1,p2,x,y,z);
	} else { 
		/* this will convert Cartesian System II coordinates to the model coords */
		_SysIII2Mag(n,p0,p1,p2,x,y,z);
	}
				
				
	/* delete temporary variables */
	delete[] x;
	delete[] y;
	delete[] z;
}
