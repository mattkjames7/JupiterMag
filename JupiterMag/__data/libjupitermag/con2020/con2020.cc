#include "con2020.h"

Con2020::Con2020() {
	/* set all model parameters to their default values */
	mui_ = 139.6;
	irho_ = 16.7;
	r0_ = 7.8;
	r0sq_ = r0_*r0_;
	r1_ = 51.4;
	r1sq_ = r1_*r1_;
	d_ = 3.6;
	xt_ = 9.3;
	xp_ = 155.8;
	strcpy(eqtype_,"hybrid");
	Edwards_ = true;
	ErrChk_ = true;
	CartIn_ = true;
	CartOut_ = true;
	deltarho_ = 1.0;
	deltaz_ = 0.01;
	smooth_ = false;
	
	
	/* some other values which will only need calculating once */
	discshift_ = (xp_-180.0)*deg2rad;
	disctilt_ = xt_*deg2rad;
	cosxp_ = cos(discshift_);
	sinxp_ = sin(discshift_);
	cosxt_ = cos(disctilt_);
	sinxt_ = sin(disctilt_);

	/* initialize some things used for integration */
	_InitIntegrals();

	/* set appropriate coord conversion functions*/
	_SetIOFunctions();
	_SetModelFunctions();

}

Con2020::Con2020(double mui, double irho, double r0, double r1,
				double d, double xt, double xp, const char *eqtype,
				bool Edwards, bool ErrChk, bool CartIn, bool CartOut) {
	/* set non-boolean model parameters to their default values */
	mui_ = 139.6;
	irho_ = 16.7;
	r0_ = 7.8;
	r0sq_ = r0_*r0_;
	r1_ = 51.4;
	r1sq_ = r1_*r1_;
	d_ = 3.6;
	xt_ = 9.3;
	xp_ = 155.8;
	strcpy(eqtype_,"hybrid");
	Edwards_ = Edwards;
	ErrChk_ = ErrChk;
	CartIn_ = CartIn;
	CartOut_ = CartOut;
	deltarho_ = 1.0;
	deltaz_ = 0.01;
	smooth_ = false;
		
	/* apply custom values if they are valid */
	SetCurrentDensity(mui);
	SetRadCurrentDensity(irho);
	SetR0(r0);
	SetR1(r1);
	SetCSHalfThickness(d);
	SetCSTilt(xt);
	SetCSTiltAzimuth(xp);
	
	
	/* some other values which will only need calculating once */
	discshift_ = (xp_-180.0)*deg2rad;
	disctilt_ = xt_*deg2rad;

	/* initialize some things used for integration */
	_InitIntegrals();
	
	/* set appropriate coord conversion functions*/
	_SetIOFunctions();
	_SetModelFunctions();

}
Con2020::~Con2020() {
	_DeleteIntegrals();
}


void Con2020::_SetIOFunctions() {
	
	if (CartIn_) {
		_ConvInput = &Con2020::_SysIII2Mag;
	} else {
		_ConvInput = &Con2020::_PolSysIII2Mag;
	}
	if (CartOut_) {
		_ConvOutput = &Con2020::_BMag2SysIII;
	} else {
		_ConvOutput = &Con2020::_BMag2PolSysIII;
	}	
}

void Con2020::_SetModelFunctions() {

	/* firstly set whether we are using the Edwards or Connerney 
	 * analytical functions */
	if (Edwards_) {
		if (smooth_) {
			_LargeRho = &Con2020::_LargeRhoEdwardsSmooth;
		} else {
			_LargeRho = &Con2020::_LargeRhoEdwards;
		}
		_SmallRho = &Con2020::_SmallRhoEdwards;
	} else {
		_LargeRho = &Con2020::_LargeRhoConnerney;
		_SmallRho = &Con2020::_SmallRhoConnerney;
	}
	
	/* now we need to set which model functions we will use
	 * (analytic, intergral or hybrid) */
	if (strcmp(eqtype_,"analytic") == 0) {
		if (smooth_) {
			_Model = &Con2020::_AnalyticSmooth;
		} else {
			_Model = &Con2020::_Analytic;
		}
	} else if (strcmp(eqtype_,"integral") == 0) {
		_Model = &Con2020::_Integral;
	} else if (strcmp(eqtype_,"hybrid") == 0) {
		_Model = &Con2020::_Hybrid;
	} else {
		printf("What's going on here then?\n");
	}
	 
}

void Con2020::_SysIII2Mag(int n, double *x0, double *y0, double *z0,
						double *x1, double *y1, double *z1, double *rho, double *absz,
						double *cost, double *sint, double *cosp, double *sinp) {
	
	
	/* some temporary variables which get used more than once */
	double xt, theta, phi, r, rho0, rho0_sq;

	
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
		xt = rho0*(cosp[i]*cosxp_ + sinp[i]*sinxp_);
		y1[i] = rho0*(sinp[i]*cosxp_ - cosp[i]*sinxp_);

		/*align with the current sheet */
		x1[i] = xt*cosxt_ + z0[i]*sinxt_;
		z1[i] = z0[i]*cosxt_ - xt*sinxt_;
		
		/* calculate rho and abs(z) */
		rho[i] = sqrt(x1[i]*x1[i] + y1[i]*y1[i]);
		absz[i] = fabs(z1[i]);
		
	}
	
					
}

void Con2020::_PolSysIII2Mag(int n, double *r, double *t, double *p,
						double *x1, double *y1, double *z1, double *rho, double *absz,
						double *cost, double *sint, double *cosp, double *sinp) {
	
	
	/* some temporary variables which get used more than once */
	double  x, z;
	
	int i;
	for (i=0;i<n;i++) {
		/* temporary variables */
		sint[i] = sin(t[i]);
		cost[i] = cos(t[i]);
		sinp[i] = sin(p[i]);
		cosp[i] = cos(p[i]);
		
		/*faster conversion code from con2020 */
		x = r[i]*sint[i]*(cosp[i]*cosxp_ + sinp[i]*sinxp_);
		y1[i] = r[i]*sint[i]*(sinp[i]*cosxp_ - cosp[i]*sinxp_);
		z = r[i]*cost[i];
		
		/*newly rotated coords */
		x1[i] = x*cosxt_ + z*sinxt_;
		z1[i] = z*cosxt_ - x*sinxt_;
		
		rho[i] = sqrt(x1[i]*x1[i] + y1[i]*y1[i]);
		absz[i] = fabs(z1[i]);
		
	}
	
					
}

void Con2020::_BMag2SysIII(int n, double *x1, double *y1, double *rho1,
				double *cost, double *sint, double *cosp, double *sinp,
				double *Brho1, double *Bphi1, double *Bz1,
				double *Bx0, double *By0, double *Bz0) {

	double cosp1, sinp1, Bx, Bx1, By1;
	
	int i;
	for (i=0;i<n;i++) {
		/* rotating longitude */
		cosp1 = x1[i]/rho1[i];
		sinp1 = y1[i]/rho1[i];
		Bx1 = Brho1[i]*cosp1 - Bphi1[i]*sinp1;
		By1 = Brho1[i]*sinp1 + Bphi1[i]*cosp1;
		
		/* rotate back by dipole tilt */
		Bx = Bx1*cosxt_ - Bz1[i]*sinxt_;
		Bz0[i] = Bx1*sinxt_ + Bz1[i]*cosxt_;
		
		/* shift to System III */
		Bx0[i] = Bx*cosxp_ - By1*sinxp_;
		By0[i] = By1*cosxp_ + Bx*sinxp_;
		
	}
	
}

void Con2020::_BMag2PolSysIII(int n, double *x1, double *y1, double *rho1,
				double *cost, double *sint, double *cosp, double *sinp,
				double *Brho1, double *Bphi1, double *Bz1,
				double *Br0, double *Bt0, double *Bp0) {

	double cosp1, sinp1, Bx, Bz, Bx1, By1, Bx2, By2;


	int i;
	for (i=0;i<n;i++) {
		
		/* rotating longitude */
		cosp1 = x1[i]/rho1[i];
		sinp1 = y1[i]/rho1[i];
		Bx1 = Brho1[i]*cosp1 - Bphi1[i]*sinp1;
		By1 = Brho1[i]*sinp1 + Bphi1[i]*cosp1;
		
		/* rotate back by dipole tilt */
		Bx = Bx1*cosxt_ - Bz1[i]*sinxt_;
		Bz = Bx1*sinxt_ + Bz1[i]*cosxt_;
		
		/* shift to System III */
		Bx2 = Bx*cosxp_ - By1*sinxp_;
		By2 = By1*cosxp_ + Bx*sinxp_;
		
		/* convert to polar */
		Br0[i] =  Bx2*sint[i]*cosp[i] + By2*sint[i]*sinp[i] + Bz*cost[i];
		Bt0[i] =  Bx2*cost[i]*cosp[i] + By2*cost[i]*sinp[i] - Bz*sint[i];
		Bp0[i] = -Bx2*sinp[i] + By2*cosp[i];
		
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

void Con2020::_AzimuthalField(double rho, double absz, double z, double *Bphi) {
	
	Bphi[0] = 2.7975*irho_/rho;
		
	if (absz < d_) {
		Bphi[0] = Bphi[0]*absz/d_;
	}
		
	if (z > 0.0) {
		Bphi[0] = -Bphi[0];
	}
}


void Con2020::Field(double p0, double p1, double p2,
					double *B0, double *B1, double *B2) {

	/* create a bunch of empty variables */
	double x, y, z, rho, absz;
	double cost, sint, cosp, sinp;	
	double Brho, Bphi, Bz;
	
	/* convert the input coordinates */
	(this->*_ConvInput)(1,&p0,&p1,&p2,&x,&y,&z,&rho,&absz,&cost,&sint,&cosp,&sinp);
	
	/* calculate the mode field */
	(this->*_Model)(rho,absz,z,&Brho,&Bphi,&Bz);
	
	/*convert the output field */
	(this->*_ConvOutput)(1,&x,&y,&rho,&cost,&sint,&cosp,&sinp,&Brho,&Bphi,&Bz,B0,B1,B2);
}


void Con2020::Field(int n, double *p0, double *p1, double *p2, 
			double *B0, double *B1, double *B2) {

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

	/* convert input coordinates to the coordinate system used by the model */
	(this->*_ConvInput)(n,p0,p1,p2,x,y,z,rho,absz,cost,sint,cosp,sinp);

	/* call the model */
	for (i=0;i<n;i++) {
		(this->*_Model)(rho[i],absz[i],z[i],&Brho[i],&Bphi[i],&Bz[i]);
	}

	/* convert B field back to appropriate coordinate system */
	(this->*_ConvOutput)(n,x,y,rho,cost,sint,cosp,sinp,Brho,Bphi,Bz,B0,B1,B2);
				
	/* delete temporary variables */
	delete[] x;
	delete[] y;
	delete[] z;
	delete[] absz;
	delete[] rho;
	delete[] sint;
	delete[] sinp;
	delete[] cost;
	delete[] cosp;
	delete[] Brho;
	delete[] Bphi;
	delete[] Bz;


}

void Con2020::_Analytic(double rho, double absz, double z, 
						double *Brho, double *Bphi, double *Bz) {
	
	/* calculate the inner contribution to Brho and Bz */
	_AnalyticInner(rho,z,Brho,Bz);
	
	/* also the azimuthal field */
	_AzimuthalField(rho,absz,z,Bphi);
	
	/* we need to calculate the outer edge contribution */
	double oBrho, oBz;
	_AnalyticOuter(rho,z,&oBrho,&oBz);
	
	/*subtract it */
	Brho[0] -= oBrho;
	Bz[0] -= oBz;
	
}

void Con2020::_AnalyticSmooth(double rho, double absz, double z, 
						double *Brho, double *Bphi, double *Bz) {
	
	/* calculate the inner contribution to Brho and Bz */
	_AnalyticInnerSmooth(rho,z,Brho,Bz);
	
	/* also the azimuthal field */
	_AzimuthalField(rho,absz,z,Bphi);
	
	/* we need to calculate the outer edge contribution */
	double oBrho, oBz;
	_AnalyticOuterSmooth(rho,z,&oBrho,&oBz);
	
	/*subtract it */
	Brho[0] -= oBrho;
	Bz[0] -= oBz;
	
}

void Con2020::_AnalyticInner(	double rho, double z, 
								double *Brho, double *Bz) {
	
	/* define a few required variables */
	double zpd = z + d_;
	double zmd = z - d_;

	if (rho >= r0_) {
		(this->*_LargeRho)(rho,z,zmd,zpd,r0sq_,Brho,Bz);
	} else {
		(this->*_SmallRho)(rho,z,zmd,zpd,r0sq_,Brho,Bz);
	}
									
}

void Con2020::_AnalyticInnerSmooth(	double rho, double z, 
									double *Brho, double *Bz) {
	
	/* define a few required variables */
	double zpd = z + d_;
	double zmd = z - d_;
	
	/* define some temporary variables */
	double Brho0, Brho1, Bz0, Bz1;
	double tanhrho, C0, C1; 
	
	/* calculate small and large rho approximations
	 * NOTE: Add smoothed functions for small and large rho */
	(this->*_LargeRho)(rho,z,zmd,zpd,r0sq_,&Brho1,&Bz1);
	(this->*_SmallRho)(rho,z,zmd,zpd,r0sq_,&Brho0,&Bz0);
	
	/* calculate the tanh smoothing parameters */
	tanhrho = tanh((rho-r0_)/deltarho_);
	C0 = (1-tanhrho)/2.0;
	C1 = (1+tanhrho)/2.0;
	
	/* splice together as suggested by Stan */
	*Brho = Brho0*C0 + Brho1*C1;
	*Bz = Bz0*C0 + Bz1*C1;
	
									
}

void Con2020::_AnalyticOuter(	double rho, double z, 
								double *Brho, double *Bz) {
	
	/* define a few required variables */
	double zpd = z + d_;
	double zmd = z - d_;

	if (rho >= r1_) {
		(this->*_LargeRho)(rho,z,zmd,zpd,r1sq_,Brho,Bz);
	} else {
		(this->*_SmallRho)(rho,z,zmd,zpd,r1sq_,Brho,Bz);
	}
									
}


void Con2020::_AnalyticOuterSmooth(	double rho, double z, 
									double *Brho, double *Bz) {
	
	/* define a few required variables */
	double zpd = z + d_;
	double zmd = z - d_;
	
	/* define some temporary variables */
	double Brho0, Brho1, Bz0, Bz1;
	double tanhrho, C0, C1; 
	
	/* calculate small and large rho approximations
	 * NOTE: Add smoothed functions for small and large rho */
	(this->*_LargeRho)(rho,z,zmd,zpd,r1sq_,&Brho1,&Bz1);
	(this->*_SmallRho)(rho,z,zmd,zpd,r1sq_,&Brho0,&Bz0);
	
	/* calculate the tanh smoothing parameters */
	tanhrho = tanh((rho-r1_)/deltarho_);
	C0 = (1-tanhrho)/2.0;
	C1 = (1+tanhrho)/2.0;
	
	/* splice together as suggested by Stan */
	*Brho = Brho0*C0 + Brho1*C1;
	*Bz = Bz0*C0 + Bz1*C1;
	
									
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

void Con2020::_LargeRhoEdwardsSmooth(double rho, double z, double zmd,
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
	double terma2 = (2.0/rho)*smoothd(z,deltaz_,d_);
	Brho[0] = mui_*(terma0 + terma1 + terma2);

	/* equation 13b */
	double termb0 = log((zpd + f2)/(zmd + f1));
	double termb1 = (a2/4.0)*(zpd/f2cubed - zmd/f1cubed);
	Bz[0] = mui_*(termb0 + termb1);
}


void Con2020::_SmallRhoEdwards(double rho, double z, double zmd, double zpd, 
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
		rnbes_[zcase] = (int) (rlmx_array_[zcase]/dlambda_brho_) - 1;
		znbes_[zcase] = (int) (zlmx_array_[zcase]/dlambda_bz_) - 1;
		
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
	}
	
	_RecalcIntegrals();
}

void Con2020::_RecalcIntegrals(){
	
	int zcase, i;
	double ld;
	for (zcase=0;zcase<6;zcase++) {	
		/* initialize lambda for z and rho cases */
		for (i=0;i<rnbes_[zcase];i++) {
			rlambda_[zcase][i] = (i+1)*dlambda_brho_;
		}
		for (i=0;i<znbes_[zcase];i++) {
			zlambda_[zcase][i] = (i+1)*dlambda_bz_;
		}
	
		/* calculate j0(r0*lambda) */
		j0(rnbes_[zcase],&rlambda_[zcase][0],r0_,&rj0_lambda_r0_[zcase][0]);
		j0(znbes_[zcase],&zlambda_[zcase][0],r0_,&zj0_lambda_r0_[zcase][0]);
		
		/* equations */
		for (i=0;i<rnbes_[zcase];i++) {
			ld = d_*rlambda_[zcase][i];
			Eq14_[zcase][i] = rj0_lambda_r0_[zcase][i]*sinh(ld)/rlambda_[zcase][i];
			Eq17_[zcase][i] = rj0_lambda_r0_[zcase][i]*exp(-ld)/rlambda_[zcase][i];
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

void Con2020::_Integral(double rho, double absz, double z, 
						double *Brho, double *Bphi, double *Bz) {
	
	/* calculate the inner contribution to Brho and Bz */
	_IntegralInner(rho,absz,z,Brho,Bz);
	
	/* also the azimuthal field */
	_AzimuthalField(rho,absz,z,Bphi);

	/* we need to calculate the outer edge contribution */
	double oBrho, oBz;
	_AnalyticOuter(rho,z,&oBrho,&oBz);

	/*subtract it */
	Brho[0] -= oBrho;
	Bz[0] -= oBz;
	
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

void Con2020::_IntegralCheck(double absz, int *chind) {

	int i;
	double check1;
	bool check2;
	
	
	check1 = fabs(absz - d_);
	check2 = (absz < d_*1.1);

	if (check1 >= 0.7) {
		chind[0] = 1;
	} else if (check1 < 0.1) {
		chind[0] = 5;
	} else {
		chind[0] = 3;
	}
	chind[0] -= ((int) check2);
	
	
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

void Con2020::_IntegralInner( 	double rho, double absz, double z,
								double *Brho, double *Bz) {
	//printf("\n_IntegralInner\n");
	int chind;
	/* check which set of integral parameters we need to use*/
	_IntegralCheck(absz,&chind);
	
	if (absz > d_) {
		/* in this case we want to use equations 14 and 15 outside
		 * of the finite current sheet */
		_IntegrateEq14(chind,rho,z,absz,Brho);
		_IntegrateEq15(chind,rho,absz,Bz);
	} else { 
		/* here we are inside the current sheet and should integrate
		 * equations 17 and 18 */
		_IntegrateEq17(chind,rho,z,Brho);
		_IntegrateEq18(chind,rho,z,Bz);			
	}	
	
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
	double *func = new double[n]; //these arrays can be pre-allocated for overwriting
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
	
					
void Con2020::_Hybrid(double rho, double absz, double z, 
						double *Brho, double *Bphi, double *Bz) {


	/* calculate the inner contribution to Brho and Bz */
	if ((absz <= 1.5*d_) && (fabs(rho - r0_) <= 2.0)) {
		/* use integration */
		_IntegralInner(rho,absz,z,Brho,Bz);
	} else {
		/* use analytical */
		_AnalyticInner(rho,z,Brho,Bz);
	}

	/* also the azimuthal field */
	_AzimuthalField(rho,absz,z,Bphi);

	/* we need to calculate the outer edge contribution */
	double oBrho, oBz;
	_AnalyticOuter(rho,z,&oBrho,&oBz);

	/*subtract it */
	Brho[0] -= oBrho;
	Bz[0] -= oBz;

}


void Con2020::SetEdwardsEqs(bool Edwards) {
	Edwards_ = Edwards;
	
	/* reset function pointers */
	_SetModelFunctions();
}

bool Con2020::GetEdwardsEqs() {
	return Edwards_;
}

void Con2020::SetEqType(const char *eqtype) {
	
	if (	(strcmp(eqtype,"analytic") == 0) ||
			(strcmp(eqtype,"integral") == 0) ||
			(strcmp(eqtype,"hybrid") == 0)) {
		/* this is a valid string - update it */
		strcpy(eqtype_,eqtype);
		
		_SetModelFunctions();
	} else {
		printf("eqtype '%s' not recognised - ignoring\n",eqtype);
	}

}

void Con2020::SetSmooth(bool smooth) {
	smooth_ = smooth;
	_SetModelFunctions();
}

bool Con2020::GetSmooth() {
	
	return smooth_;
}
	

void Con2020::GetEqType(char *eqtype) {
	strcpy(eqtype,eqtype_);
}
	
void Con2020::SetCurrentDensity(double mui) {
	
	if (isfinite(mui)) {
		/* good value (hopefully) */
		mui_ = mui;
	} else {
		printf("Non-finite value - ignoring\n");
	}
}

double Con2020::GetCurrentDensity() {
	return mui_;
}
	
void Con2020::SetRadCurrentDensity(double irho) {
	if (isfinite(irho)) {
		/* good value (hopefully) */
		irho_ = irho;
	} else {
		printf("Non-finite value - ignoring\n");
	}	
}

double Con2020::GetRadCurrentDensity() {
	return irho_;
}

void Con2020::SetR0(double r0) {
	if (isfinite(r0) && (r0 >= 0.0)) {
		/* good value (hopefully) */
		r0_ = r0;
		r0sq_ = r0_*r0_;
	} else if (!isfinite(r0)) {
		printf("Non-finite value - ignoring\n");
	} else {
		printf("r0 must have a positive value\n");
	}
}

double Con2020::GetR0() {
	return r0_;
}

void Con2020::SetR1(double r1) {
	if (isfinite(r1) && (r1 >= 0.0)) {
		/* good value (hopefully) */
		r1_ = r1;
		r1sq_ = r1_*r1_;
	} else if (!isfinite(r1)) {
		printf("Non-finite value - ignoring\n");
	} else {
		printf("r1 must have a positive value\n");
	}
}

double Con2020::GetR1() {
	return r1_;
}

void Con2020::SetCSHalfThickness(double d) {

	if (isfinite(d) && (d >= 0.0)) {
		/* good value (hopefully) */
		d_ = d;
	} else if (!isfinite(d)) {
		printf("Non-finite value - ignoring\n");
	} else {
		printf("d must have a positive value\n");
	}
}	

double Con2020::GetCSHalfThickness() {
	return d_;
}
	
void Con2020::SetCSTilt(double xt) {
	
	if (isfinite(xt)) {
		/* good value (hopefully) */
		xt_ = xt;
		disctilt_ = xt_*deg2rad;
		cosxt_ = cos(disctilt_);
		sinxt_ = sin(disctilt_);
	} else {
		printf("Non-finite value - ignoring\n");
	}
}

double Con2020::GetCSTilt() {
	return xt_;
}

void Con2020::SetCSTiltAzimuth(double xp) {
	if (isfinite(xp)) {
		/* good value (hopefully) */
		xp_ = xp;
		discshift_ = (xp_ - 180.0)*deg2rad;
		cosxp_ = cos(discshift_);
		sinxp_ = sin(discshift_);
	} else {
		printf("Non-finite value - ignoring\n");
	}
}		
		
double Con2020::GetCSTiltAzimuth() {
	return xp_;
}
		
void Con2020::SetErrCheck(bool ErrChk) {
	ErrChk_ = ErrChk;
}

bool Con2020::GetErrCheck() {
	return ErrChk_;
}

void Con2020::SetCartIn(bool CartIn) {
	CartIn_ = CartIn;
	_SetIOFunctions();
}

bool Con2020::GetCartIn() {
	return CartIn_;
}

void Con2020::SetCartOut(bool CartOut) {
	CartOut_ = CartOut;
	_SetIOFunctions();
}

bool Con2020::GetCartOut() {
	return CartOut_;
}
