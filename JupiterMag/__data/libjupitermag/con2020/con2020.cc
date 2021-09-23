#include "con2020.h"

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
	discshift_ = xp_*deg2rad;
	disctilt_ = xt_*deg2rad;

	/* initialize some things used for integration */
	_InitIntegrals();

}
Con2020::~Con2020() {
	_DeleteIntegrals();
}

void Con2020::_SysIII2Mag(int n, double *x0, double *y0, double *z0,
						double *x1, double *y1, double *z1, double *rho, double *absz,
						double *sint, double *cost, double *sinp, double *cosp) {
	
	
	/* some temporary variables which get used more than once */
	double sincst, coscst, coscss, sincss, xt, theta, phi, r, rho0, rho0_sq;
	sincss = sin(discshift_);
	coscss = cos(discshift_);
	sincst = sin(disctilt_);
	coscst = cos(disctilt_);

	
	int i;
	for (i=0;i<n;i++) {
		/*calculate some angles and stuff*/
		rho0_sq = x0[i]*x0[i] + y0[i]*y0[i];
		rho0 = sqrt(rho0_sq);
		r = sqrt(rho0_sq + z0[i]*z0[i]);
		
		cost[i] = z0[i]/r;
		sint[i] = rho0/r;
		sinp[i] = y0[i]/rho0;
		cosp[i] = x0[i]/rho0;
		
		/*rotate about z0 for the correct longitude */
		xt = rho0*(cosp[i]*coscss + sinp[i]*sincss);
		y1[i] = rho0*(sinp[i]*coscss - cosp[i]*sincss);

		/*align with the current sheet */
		x1[i] = xt*coscst + z0[i]*sincst;
		z1[i] = z0[i]*coscst - xt*sincst;
		
		/* calculate rho and abs(z) */
		rho[i] = sqrt(x1[i]*x1[i] + y1[i]*y1[i]);
		absz[i] = fabs(z1[i]);
		
	}
	
					
}

void Con2020::_PolSysIII2Mag(int n, double *r, double *t, double *p,
						double *x1, double *y1, double *z1, double *rho, double *absz,
						double *sint, double *cost, double *sinp, double *cosp) {
	
	
	/* some temporary variables which get used more than once */
	double sindt, cosdt, cospds, sinpds, rsint, rcost;

	sindt = sin(disctilt_);
	cosdt = cos(disctilt_);
	
	
	int i;
	for (i=0;i<n;i++) {
		/* temporary variables */
		sint[i] = sin(t[i]);
		cost[i] = cos(t[i]);
		sinp[i] = sin(p[i]);
		cosp[i] = cos(p[i]);
		rsint = r[i]*sint[i];
		rcost = r[i]*cost[i];
		sinpds = sin(p[i] - discshift_);
		cospds = cos(p[i] - discshift_);		
		
		/*newly rotated coords */
		x1[i] = rsint*cospds*cosdt + rcost*sindt;
		y1[i] = rsint*sinpds;
		z1[i] = rcost*cosdt - rsint*cospds*sindt;
		
		rho[i] = sqrt(x1[i]*x1[i] + y1[i]*y1[i]);
		absz[i] = fabs(z1[i]);
		
	}
	
					
}

void Con2020::_BMag2SysIII(int n, double *x, double *y, double *rho,
				double *Brho, double *Bphi, double *Bz,
				double *Bxo, double *Byo, double *Bzo) {

	double sindt, cosdt, sinds, cosds, cospl, sinpl;
	double Bx1, By1, Bx2;

	sindt = sin(disctilt_);
	cosdt = cos(disctilt_);
	sinds = sin(discshift_);
	cosds = cos(discshift_);
	
	int i;
	for (i=0;i<n;i++) {
		
		/* rotating longitude */
		cospl = x[i]/rho[i];
		sinpl = y[i]/rho[i];
		Bx1 = Brho[i]*cospl - Bphi[i]*sinpl;
		By1 = Brho[i]*sinpl + Bphi[i]*cospl;
		
		/* rotate back by dipole tilt */
		Bx2 = Bx1*cosdt - Bz[i]*sindt;
		Bzo[i] = Bx1*sindt + Bz[i]*cosdt;
		
		/* shift to System III */
		Bxo[i] = Bx2*cosds - By1*sinds;
		Byo[i] = By1*cosds + Bx2*sinds;
		
	}
	
}

void Con2020::_BMag2PolSysIII(int n, double *x, double *y, double *rho,
				double *sint, double *cost, double *sinp, double *cosp,
				double *Brho, double *Bphi, double *Bz,
				double *Br, double *Bt, double *Bp) {

	double sindt, cosdt, sinds, cosds, cospl, sinpl;
	double Bx1, By1, Bx2, Bz2, Bx3, By3;

	sindt = sin(disctilt_);
	cosdt = cos(disctilt_);
	sinds = sin(discshift_);
	cosds = cos(discshift_);
	
	int i;
	for (i=0;i<n;i++) {
		
		/* rotating longitude */
		cospl = x[i]/rho[i];
		sinpl = y[i]/rho[i];
		Bx1 = Brho[i]*cospl - Bphi[i]*sinpl;
		By1 = Brho[i]*sinpl + Bphi[i]*cospl;
		
		/* rotate back by dipole tilt */
		Bx2 = Bx1*cosdt - Bz[i]*sindt;
		Bz2 = Bx1*sindt + Bz[i]*cosdt;
		
		/* shift to System III */
		Bx3 = Bx2*cosds - By1*sinds;
		By3 = By1*cosds + Bx2*sinds;
		
		/* convert to polar */
		Br[i] = Bx3*sint[i]*cosp[i] + By3*sint[i]*sinp[i] + Bz2*cost[i];
		Bt[i] = Bx3*cost[i]*cosp[i] + By3*cost[i]*sint[i] + Bz2*sint[i];
		Bp[i] = -Bx3*sinp[i] + By3*cosp[i];
		
	}
	
}

void Con2020::_AzimuthalField(int n, double *rho, double *z, double *absz, double *Bphi) {
	
	int i;
	for (i=0;i<n;i++) {
		Bphi[i] = 2.7975*irho_/rho[i];
		
		if (absz[i] < d_) {
			Bphi[i] = Bphi[i]*absz[i]/d_;
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
	double *absz = new double[n];
	double *sint = new double[n];
	double *sinp = new double[n];
	double *cost = new double[n];
	double *cosp = new double[n];
	double *rho = new double[n];
	double *Brho = new double[n];
	double *Bphi = new double[n];
	double *Bz = new double[n];
	double *Brhoo = new double[n];
	double *Bzo = new double[n];

	/* convert input coordinates to the coordinate system used by the model */
	if (PolIn) {
		/* this will convert spherical polar to the Cartesian model coordinates */
		_PolSysIII2Mag(n,p0,p1,p2,x,y,z,rho,absz,sint,cost,sinp,cosp);
	} else { 
		/* this will convert Cartesian System III coordinates to the model coords */
		_SysIII2Mag(n,p0,p1,p2,x,y,z,rho,absz,sint,cost,sinp,cosp);
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
			if ((fabs(z[i]) <= d15) && (fabs(rho[i]-r0_) <= 2.0)) {
				indint[nint] = i;
				nint++;
			} else { 
				indana[nana] = i;
				nana++;
			}
		}
	}
	
	/* more temporary position arrays */
	double *rhoi, *zi, *abszi, *Brhoi, *Bzi;	
	
	/* call integral function if needed */
	if (nint > 0) {
		/* allocate temporary variables */
		rhoi = new double[nint];
		zi = new double[nint];
		abszi = new double[nint];
		Brhoi = new double[nint];
		Bzi = new double[nint];
		
		/* move positions into new array */
		for (i=0;i<nint;i++) {
			rhoi[i] = rho[indint[i]];
			zi[i] = z[indint[i]];
			abszi[i] = absz[indint[i]];
		}
		
		/* call the integral model */
		_SolveIntegral(nint,rhoi,zi,abszi,Brhoi,Bzi);
		
		/* move into the new array */
		for (i=0;i<nint;i++) {
			Brho[indint[i]] = Brhoi[i];
			Bz[indint[i]] = Bzi[i];
		}
		
		/* deallocate temporary vars */
		delete[] rhoi;
		delete[] zi;
		delete[] abszi;
		delete[] Brhoi;
		delete[] Bzi;
	}
	
	/* call analytical fucntion if needed */
	if (nana > 0) {
		/* allocate temporary variables */
		rhoi = new double[nana];
		zi = new double[nana];
		Brhoi = new double[nana];
		Bzi = new double[nana];
		
		/* move positions into the new array */
		for (i=0;i<nana;i++) {
			rhoi[i] = rho[indana[i]];
			zi[i] = z[indana[i]];
		}
		
		/* call the analytical model */
		_SolveAnalytic(nana,rhoi,zi,r0_,Brhoi,Bzi);
		
		/* move into the new array */
		for (i=0;i<nana;i++) {
			Brho[indana[i]] = Brhoi[i];
			Bz[indana[i]] = Bzi[i];
		}
		
		/* deallocate temporary vars */
		delete[] rhoi;
		delete[] zi;
		delete[] Brhoi;
		delete[] Bzi;
	}	

	/* calculate the outer edge of the current sheet */
	_SolveAnalytic(n,rho,z,r1_,Brhoo,Bzo);	
	
	/* subtract outer contribution from inner one */
	for (i=0;i<n;i++) {
		Brho[i] = Brho[i] - Brhoo[i];
		Bz[i] = Bz[i] - Bzo[i];
	}
	
	/* calculate azimuthal field */	
	_AzimuthalField(n,rho,z,absz,Bphi);
	
		
	/* convert B field back to appropriate coordinate system */
	if (PolOut) {
		/* convert to the polar coordinate system (it's actually the
		 * Cartesian components of the field oriented such that there
		 * is a radial component, phi component and theta component - as
		 * opposed to a magnitude and two angles) */
		 _BMag2PolSysIII(n,x,y,rho,sint,cost,sinp,cosp,Brho,Bphi,Bz,B0,B1,B2);
	} else {
		/* convert to System III Cartesian */
		_BMag2SysIII(n,x,y,rho,Brho,Bphi,Bz,B0,B1,B2);
	}
				
				
	/* delete temporary variables */
	delete[] x;
	delete[] y;
	delete[] z;
	delete[] absz;
	delete[] rho;
	delete[] indint;
	delete[] indana;
	delete[] sint;
	delete[] sinp;
	delete[] cost;
	delete[] cosp;
	delete[] Brho;
	delete[] Bphi;
	delete[] Bz;
	delete[] Brhoo;
	delete[] Bzo;

	
}


void Con2020::_SolveAnalytic(int n, double *rho, double *z, double a, 
				double *Brho, double *Bz) {
	
	int i;
	double zpd, zmd, a2;
	a2 = a*a;
	if (Edwards_) {
		/* solve the Edwards et al equations */
		for (i=0;i<n;i++) {
			zpd = z[i] + d_;
			zmd = z[i] - d_;
			if (rho[i] >= a) {
				/* large rho approximation */
				_LargeRhoEdwards(rho[i],z[i],zmd,zpd,a2,&Brho[i],&Bz[i]);
			} else {
				/* small rho approximation */
				_SmallRhoEdwards(rho[i],zmd,zpd,a2,&Brho[i],&Bz[i]);
			}
		}
	} else { 
		/* Solve the Connerney et al equations */
		for (i=0;i<n;i++) {
			zpd = z[i] + d_;
			zmd = z[i] - d_;
			if (rho[i] >= a) {
				/* large rho approximation */
				_LargeRhoConnerney(rho[i],z[i],zmd,zpd,a2,&Brho[i],&Bz[i]);
			} else {
				/* small rho approximation */
				_SmallRhoConnerney(rho[i],z[i],zmd,zpd,a2,&Brho[i],&Bz[i]);
			}
		}		
	}
}

void Con2020::_LargeRhoConnerney(double rho, double z, double zmd, 
					double zpd, double a2, double *Brho, double *Bz) {

	
	/* some common variables */
	double zmd2 = zmd*zmd;
	double zpd2 = zpd*zpd;
	double rho2 = rho*rho;
	double f1 = sqrt(zmd2 + rho2);
	double f2 = sqrt(zpd2 + rho2);
	double f1cubed = f1*f1*f1;
	double f2cubed = f2*f2*f2;
	
	/* Equation A7 */
	double termr0 = (1.0/rho)*(f1 -f2 + 2*clip(z,-d_,d_));
	double termr1 = (a2*rho/4.0)*(1/f1cubed - 1/f2cubed);
	Brho[0] = mui_*(termr0 - termr1);
	
	/* Equation A8 */
	double termz0 = 2*d_/sqrt(z*z + rho*rho);
	double termz1 = (a2/4.0)*(zmd/f1cubed - zpd/f2cubed);
	Bz[0] = mui_*(termz0 - termz1);
		
}

void Con2020::_SmallRhoConnerney(double rho, double z, double zmd, double zpd, 
								double a2, double *Brho, double *Bz) {

	
	double zmd2 = zmd*zmd;
	double zpd2 = zpd*zpd;
	double f1 = sqrt(zmd2 + a2);
	double f2 = sqrt(zpd2 + a2);
	double f1cubed = f1*f1*f1;
	double f2cubed = f2*f2*f2;

	/* Equation A1 */
	Brho[0] = mui_*(rho/2.0)*(1.0/f1 - 1.0/f2);
	
	/* Equation A2 */
	Bz[0] = mui_*(2*d_*(1/sqrt(z*z + a2)) - ((rho*rho)/4)*((zmd/f1cubed) - (zpd/f2cubed)));

}


void Con2020::_LargeRhoEdwards(double rho, double z, double zmd,
					double zpd, double a2,double *Brho, double *Bz) {
	
	
	/* some common variables */
	double zmd2 = zmd*zmd;
	double zpd2 = zpd*zpd;
	double rho2 = rho*rho;
	double f1 = sqrt(zmd2 + rho2);
	double f2 = sqrt(zpd2 + rho2);
	double f1cubed = f1*f1*f1;
	double f2cubed = f2*f2*f2;
	
	/* equation 13a */
	double terma0 = (1.0/rho)*(f1 - f2);
	double terma1 = (rho*a2/4)*(1.0/f2cubed - 1.0/f1cubed);
	double terma2 = (2.0/rho)*clip(z,-d_,d_);
	Brho[0] = mui_*(terma0 + terma1 + terma2);
	
	/* equation 13b */
	double termb0 = log((zpd + f2)/(zmd + f1));
	double termb1 = (a2/4.0)*(zpd/f2cubed - zmd/f1cubed);
	Bz[0] = mui_*(termb0 + termb1);
}


void Con2020::_SmallRhoEdwards(double rho, double zmd, double zpd, 
					double a2, double *Brho, double *Bz) { 
	
	double zmd2 = zmd*zmd;
	double zpd2 = zpd*zpd;
	double f1 = sqrt(zmd2 + a2);
	double f2 = sqrt(zpd2 + a2);
	double f1cubed = f1*f1*f1;
	double f2cubed = f2*f2*f2;
	
	/* calculate some of the common terms from equations 9a and 9b */
	double rhoov2 = rho/2.0;
	double rho2ov4 = rhoov2*rhoov2;
	double rho3ov16 = rho2ov4*rhoov2/2.0;
	
	/* equation 9a */
	double f3a = f1*f1;
	double f4a = f2*f2;
	double f3 = (a2 - 2*zmd2)/(f3a*f3a*f1);
	double f4 = (a2 - 2*zpd2)/(f4a*f4a*f2);
	
	double terma0 = rhoov2*(1.0/f1 - 1.0/f2);
	double terma1 = rho3ov16*(f3 - f4);
	
	Brho[0] = mui_*(terma0 + terma1);
	
	/* equation 9b */
	double termb0 = log((zpd + f2)/(zmd + f1));
	double termb1 = rho2ov4*(zpd/f2cubed - zmd/f1cubed);
	Bz[0] = mui_*(termb0 + termb1);
}


void Con2020::_InitIntegrals() {

	/* the lambda max values for brho and bz */
	double lmx_brho, lmx_bz;
	rlmx_array_[0] = 4.0;
	rlmx_array_[1] = 4.0;
	rlmx_array_[2] = 40.0;
	rlmx_array_[3] = 40.0;
	rlmx_array_[4] = 100.0;
	rlmx_array_[5] = 100.0;
	zlmx_array_[0] = 100.0;
	zlmx_array_[1] = 20.0;
	zlmx_array_[2] = 100.0;
	zlmx_array_[3] = 20.0;
	zlmx_array_[4] = 100.0;
	zlmx_array_[5] = 20.0;

	
	/* initialize the bessel functions which do not change */
	rnbes_ = new int[6];
	znbes_ = new int[6];
	rlambda_ = new double*[6];
	zlambda_ = new double*[6];
	rj0_lambda_r0_ = new double*[6];
	zj0_lambda_r0_ = new double*[6];
	
	/*These functions do change and will need reassigning
	 * for each run of the integral function*/
	rj1_lambda_rho_ = new double*[6];
	zj0_lambda_rho_ = new double*[6];

	/*these will store bits of the Connerney et al 1981 equations which 
	 * don't change*/
	Eq14_ = new double*[6];
	Eq15_ = new double*[6];
	Eq17_ = new double*[6];
	Eq18_ = new double*[6];
	ExpLambdaD_ = new double*[6];
	

	/* initialize the second dimensions */
	int zcase, i;
	double ld;
	for (zcase=0;zcase<6;zcase++) {
		/* calculate the length of the arrays for each case */
		rnbes_[zcase] = (int) (rlmx_array_[zcase]/dlambda_brho_);
		znbes_[zcase] = (int) (zlmx_array_[zcase]/dlambda_bz_);
		
		/* allocate the second dimension */
		rlambda_[zcase] = new double[rnbes_[zcase]];
		zlambda_[zcase] = new double[znbes_[zcase]];
		rj0_lambda_r0_[zcase] = new double[rnbes_[zcase]];
		zj0_lambda_r0_[zcase] = new double[znbes_[zcase]];
		rj1_lambda_rho_[zcase] = new double[rnbes_[zcase]];
		zj0_lambda_rho_[zcase] = new double[znbes_[zcase]];
		Eq14_[zcase] = new double[rnbes_[zcase]];
		Eq15_[zcase] = new double[znbes_[zcase]];
		Eq17_[zcase] = new double[rnbes_[zcase]];
		Eq18_[zcase] = new double[znbes_[zcase]];
		ExpLambdaD_[zcase] = new double[znbes_[zcase]];
		
		/* initialize lambda for z and rho cases */
		for (i=0;i<rnbes_[zcase];i++) {
			rlambda_[zcase][i] = i*dlambda_brho_;
		}
		for (i=0;i<znbes_[zcase];i++) {
			zlambda_[zcase][i] = i*dlambda_bz_;
		}
	
		/* calculate j0(r0*lambda) */
		j0(rnbes_[zcase],&rlambda_[zcase][0],r0_,&rj0_lambda_r0_[zcase][0]);
		j0(znbes_[zcase],&zlambda_[zcase][0],r0_,&zj0_lambda_r0_[zcase][0]);
		
		/* equations */
		for (i=0;i<rnbes_[zcase];i++) {
			ld = d_*rlambda_[zcase][i];
			Eq14_[zcase][i] = rj0_lambda_r0_[zcase][i]*sinh(ld)/rlambda_[zcase][i];
			Eq17_[zcase][i] = rj0_lambda_r0_[zcase][i]*exp(-ld);
		}	
		for (i=0;i<znbes_[zcase];i++) {
			ld = d_*zlambda_[zcase][i];
			Eq15_[zcase][i] = zj0_lambda_r0_[zcase][i]*sinh(ld)/zlambda_[zcase][i];
			Eq18_[zcase][i] = zj0_lambda_r0_[zcase][i]/zlambda_[zcase][i];
			ExpLambdaD_[zcase][i] = exp(-ld);
		}
	}
	
	
		
}

void Con2020::_DeleteIntegrals() {
	
	int i;
	for (i=0;i<6;i++) {
		delete[] rlambda_[i];
		delete[] zlambda_[i];
		delete[] rj0_lambda_r0_[i];
		delete[] rj1_lambda_rho_[i];
		delete[] zj0_lambda_r0_[i];
		delete[] zj0_lambda_rho_[i];
		delete[] Eq14_[i];
		delete[] Eq15_[i];
		delete[] Eq17_[i];
		delete[] Eq18_[i];
		delete[] ExpLambdaD_[i];
	}
	
	delete[] rlambda_;
	delete[] zlambda_;
	delete[] rj0_lambda_r0_;
	delete[] rj1_lambda_rho_;
	delete[] zj0_lambda_r0_;
	delete[] zj0_lambda_rho_;
	delete[] Eq14_;
	delete[] Eq15_;
	delete[] Eq17_;
	delete[] Eq18_;
	delete[] ExpLambdaD_;
	delete[] rnbes_;
	delete[] znbes_;
}


void Con2020::_IntegralChecks(int n, double *absz, int *chind, int ncase[]) {

	int i;
	double check1;
	bool check2;
	
	for (i=0;i<6;i++) {
		ncase[i] = 0;
	}
	
	for (i=0;i<n;i++) {
		check1 = fabs(absz[i] - d_);
		check2 = (absz[i] < d_*1.1);

		if (check1 >= 0.7) {
			chind[i] = 1;
		} else if (check1 < 0.1) {
			chind[i] = 5;
		} else {
			chind[i] = 3;
		}
		chind[i] -= ((int) check2);
		ncase[chind[i]]++;
	}
	
}

void Con2020::_SolveIntegral(int n, double *rho, double *z,
					double *absz, double *Brho, double *Bz) {
	
	int i, zcase, ncase[6];
	

	/* calculate the "checks" */
	int *chind = new int[n];
	_IntegralChecks(n,absz,chind,ncase);

	/* loop through each position */
	double br, bz;
	for (i=0;i<n;i++) {
		if (absz[i] > d_) {
			/* in this case we want to use equations 14 and 15 outside
			 * of the finite current sheet */
			_IntegrateEq14(chind[i],rho[i],z[i],absz[i],&Brho[i]);
			_IntegrateEq15(chind[i],rho[i],absz[i],&Bz[i]);
			 
		} else { 
			/* here we are inside the current sheet and should integrate
			 * equations 17 and 18 */
			_IntegrateEq17(chind[i],rho[i],z[i],&Brho[i]);
			_IntegrateEq18(chind[i],rho[i],z[i],&Bz[i]);			
		}
	}
	
	delete[] chind;
	
}	
	
void Con2020::_IntegrateEq14(int zcase, double rho, double z, double absz, double *Brho) {
	
	/* create an array to integrate over */
	int n = rnbes_[zcase];
	double *func = new double[n];
	double *j1lr = new double[n];
	
	/* calculate the other bessel function */
	j1(n,&rlambda_[zcase][0],rho,j1lr);
	
	
	/* calculate the function */
	int i;
	for (i=0;i<n;i++) {
		func[i] = Eq14_[zcase][i]*j1lr[i]*exp(-rlambda_[zcase][i]*absz);
	}
	
	/* integrate it */
	Brho[0] = sgn(z)*2.0*mui_*trapc(n,dlambda_brho_,func);
	
	delete[] func;
	delete[] j1lr;
	
}
	
					
void Con2020::_IntegrateEq15(int zcase, double rho, double absz, double *Bz) {
	
	/* create an array to integrate over */
	int n = znbes_[zcase];
	double *func = new double[n];
	double *j0lr = new double[n];
	
	/* calculate the other bessel function */
	j0(n,&zlambda_[zcase][0],rho,j0lr);
	
	
	/* calculate the function */
	int i;
	for (i=0;i<n;i++) {
		func[i] = Eq15_[zcase][i]*j0lr[i]*exp(-zlambda_[zcase][i]*absz);
	}
	
	/* integrate it */
	Bz[0] = 2.0*mui_*trapc(n,dlambda_bz_,func);

	delete[] func;
	delete[] j0lr;
	
}
	
					
	
void Con2020::_IntegrateEq17(int zcase, double rho, double z, double *Brho) {
	
	/* create an array to integrate over */
	int n = rnbes_[zcase];
	double *func = new double[n];
	double *j1lr = new double[n];
	
	/* calculate the other bessel function */
	j1(n,&rlambda_[zcase][0],rho,j1lr);
	
	
	/* calculate the function */
	int i;
	for (i=0;i<n;i++) {
		func[i] = Eq17_[zcase][i]*j1lr[i]*sinh(rlambda_[zcase][i]*z);
	}
	
	/* integrate it */
	Brho[0] = 2.0*mui_*trapc(n,dlambda_brho_,func);

	delete[] func;
	delete[] j1lr;
	
}
	
					
void Con2020::_IntegrateEq18(int zcase, double rho, double z, double *Bz) {
	
	/* create an array to integrate over */
	int n = znbes_[zcase];
	double *func = new double[n];
	double *j0lr = new double[n];
	
	/* calculate the other bessel function */
	j0(n,&zlambda_[zcase][0],rho,j0lr);
	
	
	/* calculate the function */
	int i;
	for (i=0;i<n;i++) {
		func[i] = Eq18_[zcase][i]*j0lr[i]*(1.0 - ExpLambdaD_[zcase][i]*cosh(zlambda_[zcase][i]*z));
	}
	
	/* integrate it */
	Bz[0] = 2.0*mui_*trapc(n,dlambda_bz_,func);

	delete[] func;
	delete[] j0lr;
	
}
	
					

