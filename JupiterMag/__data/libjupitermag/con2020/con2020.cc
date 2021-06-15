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
