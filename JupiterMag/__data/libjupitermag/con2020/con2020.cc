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

void Con2020::AzimuthalField(int n, double *rho, double *z, double *zabs, double *Bphi) {
	
	int i;
	for (i=0;i<n;i++) {
		Bphi[i] = 2.7975*irho_/rho[i];
		
		if (zabs[i] < d_) {
			Bphi[i] = Bphi[i]*zabs[i]/d_;
		}
		
		if (z[i] > 0.0) {
			Bphi[i] = -Bphi[i];
		}
			
	}
}

void Con2020::Field(int n, double *p0, double *p1, double *p2, 
			double *B0, double *B1, double *B2,
			bool PolIn, bool PolOut) {

	int i;

	/* allocate arrays to store model coordinates and fields */
	double *x = new double[n];
	double *y = new double[n];
	double *z = new double[n];
	double *Brho = new double[n];
	double *Bphi = new double[n];
	double *Bz = new double[n];

	/* convert input coordinates to the coordinate system used by the model */
	if (PolIn) {
		/* this will convert spherical polar to the Cartesian model coordinates */
		_PolSysIII2Mag(n,p0,p1,p2,x,y,z);
	} else { 
		/* this will convert Cartesian System III coordinates to the model coords */
		_SysIII2Mag(n,p0,p1,p2,x,y,z);
	}
	
	/* calculate rho */
	double *rho = new double[n];
	for (i=0;i<n;i++) {
		rho[i] = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
	}
				
	/* determine which method to use for each point in the inner edge calculation */
	int nint, nana; 
	int *indint = new int[n];
	int *indana = new int[n];
	double d15 = d_*1.5;
	nint = 0;
	nana = 0;
	if (strcmp(eqtype_,"analytic") == 0) {
		/* set all to analytic */
		for (i=0;i<n;i++) {
			indana[i] = i;
			nana++;
		}
	} else if (strcmp(eqtype_,"integral") == 0) {
		/* set all to analytic */
		for (i=0;i<n;i++) {
			indint[i] = i;
			nint++;
		}		
	} else {
		/* the hybrid approach */
		for (i=0;i<n;i++) {
			if ((fabs(z[i]) <= d15) || (fabs(rho[i]-r0) <= 2.0)) {
				indint[nint] = i;
				nint++;
			} else { 
				indana[nana] = i;
				nana++;
			}
		}
	}
	
	/* more temporary position arrays */
	double *xi, *yi, *zi, *Brhoi, *Bzi;	
	
	/* call integral function if needed */
	if (nint > 0) {
		/* allocate temporary variables */
		xi = new double[nint];
		yi = new double[nint];
		zi = new double[nint];
		Brhoi = new double[nint];
		Bzi = new double[nint];
		
		/* move positions into new array */
		for (i=0;i<nint;i++) {
			xi[i] = x[indint[i]];
			yi[i] = y[indint[i]];
			zi[i] = z[indint[i]];
		}
		
		/* call the integral model */
		_SolveIntegral(nint,xi,yi,zi,Bxi,Byi,Bzi);
		
		/* move into the new array */
		for (i=0;i<nint;i++) {
			Brho[indint[i]] = Brhoi[i];
			Bz[indint[i]] = Bzi[i];
		}
		
		/* deallocate temporary vars */
		delete[] xi;
		delete[] yi;
		delete[] zi;
		delete[] Brhoi;
		delete[] Bzi;
	}
	
	/* call analytical fucntion if needed */
	if (nana > 0) {
		/* allocate temporary variables */
		xi = new double[nana];
		yi = new double[nana];
		zi = new double[nana];
		Brhoi = new double[nana];
		Bzi = new double[nana];
		
		/* move positions into the new array */
		for (i=0;i<nana;i++) {
			xi[i] = x[indana[i]];
			yi[i] = y[indana[i]];
			zi[i] = z[indana[i]];
		}
		
		/* call the analytical model */
		_SolveAnalytic(nana,xi,yi,zi,Brhoi,Bzi,r0_);
		
		/* move into the new array */
		for (i=0;i<nana;i++) {
			Brho[indana[i]] = Brhoi[i];
			Bz[indana[i]] = Bzi[i];
		}
		
		/* deallocate temporary vars */
		delete[] xi;
		delete[] yi;
		delete[] zi;
		delete[] Brhoi;
		delete[] Bzi;
	}	
	
	/* calculate the outer edge of the current sheet */
	_SolveAnalytic(n,x,y,z,Brho,Bz,r1_);	
		
	/* calculate azimuthal field */	
	_AzimuthalField(n,z,Bphi);
		
	/* convert B field back to appropriate coordinate system */
	if (PolOut) {
		/* convert to the polar coordinate system (it's actually the
		 * Cartesian components of the field oriented such that there
		 * is a radial component, phi component and theta component - as
		 * opposed to a magnitude and two angles) */
		 _BMag2Pol(n,x,y,z,rho,Brho,Bphi,Bz,B0,B1,B2);
	} else {
		/* convert to System III Cartesian */
		_BMag2SysIII(n,x,y,z,rho,Brho,Bphi,Bz,B0,B1,B2);
	}
				
				
	/* delete temporary variables */
	delete[] x;
	delete[] y;
	delete[] z;
	delete[] rho;
	delete[] indint;
	delete[] indana;
	
}
