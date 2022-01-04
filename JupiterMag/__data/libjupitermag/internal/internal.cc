#include "internal.h"

Internal vip4(&_binary_vip4coeffs_bin_start);
Internal jrm09(&_binary_jrm09coeffs_bin_start);
Internal jrm33(&_binary_jrm33coeffs_bin_start);
Internal gsfc13ev(&_binary_gsfc13evcoeffs_bin_start);
Internal gsfc15ev(&_binary_gsfc15evcoeffs_bin_start);
Internal gsfc15evs(&_binary_gsfc15evscoeffs_bin_start);
Internal isaac(&_binary_isaaccoeffs_bin_start);
Internal jpl15ev(&_binary_jpl15evcoeffs_bin_start);
Internal jpl15evs(&_binary_jpl15evscoeffs_bin_start);
Internal o4(&_binary_o4coeffs_bin_start);
Internal o6(&_binary_o6coeffs_bin_start);
Internal p11a(&_binary_p11acoeffs_bin_start);
Internal sha(&_binary_shacoeffs_bin_start);
Internal u17ev(&_binary_u17evcoeffs_bin_start);
Internal v117ev(&_binary_v117evcoeffs_bin_start);
Internal vipal(&_binary_vipalcoeffs_bin_start);
Internal vit4(&_binary_vit4coeffs_bin_start);

Internal::Internal(const char *model) {
	
	
	/* use the model string to determine which model to load */
	unsigned char *ptr;
	if (strcmp(model,"VIP4") == 0) {
		/* load the VIP 4 model */
		ptr = &_binary_vip4coeffs_bin_start;
	} else if (strcmp(model,"JRM09") == 0) {
		/* load the JRM09 model */
		ptr = &_binary_jrm09coeffs_bin_start;
	} else {
		/* default to VIP4 */
		ptr = &_binary_vip4coeffs_bin_start;
	}
	
	
	
	/* read the coeffs into the object */
	_LoadSchmidt(ptr);
	
	/* calcualte Schmidt normalized coefficient grids */
	_Schmidt();
	_CoeffGrids();
	
	/* set I/O coords */
	CartIn_ = true;
	CartOut_ = true;
	
	/* tell object that it is not a copy */
	copy = false;
	
}

Internal::Internal(unsigned char *ptr) {
	
	
	/* read the coeffs into the object */
	_LoadSchmidt(ptr);
	
	/* calcualte Schmidt normalized coefficient grids */
	_Schmidt();
	_CoeffGrids();
	
	/* set I/O coords */
	CartIn_ = true;
	CartOut_ = true;	
	
	/* tell object that it is not a copy */
	copy = false;
}

Internal::Internal(const Internal &obj) {
	
	copy = true;
	nschc_ = obj.nschc_;
	schc_ = obj.schc_;
	Snm_ = obj.Snm_;
	nmax_ = obj.nmax_;
	g_ = obj.g_;
	h_ = obj.h_;
}

Internal::~Internal() {
	if (!copy) {
		/* delete the structure containing coefficients */
		delete[] schc_;
		
		/* delete other variables containing other coefficients */
		int n;
		for (n=0;n<=nmax_;n++) {
			delete[] Snm_[n];
			delete[] g_[n];
			delete[] h_[n];
		}
		delete[] Snm_;
		delete[] g_;
		delete[] h_;
	}
}

void Internal::SetCartIn(bool CartIn) {
	CartIn_ = CartIn;
}

bool Internal::GetCartIn() {
	return CartIn_;
}

void Internal::SetCartOut(bool CartOut) {
	CartOut_ = CartOut;
}

bool Internal::GetCartOut() {
	return CartOut_;
}

void Internal::_LoadSchmidt(unsigned char *ptr){
	
	/* this is the length of each array */
	int l, i, j, p;
	
	/* read the length */
	l = ((int*) ptr)[0];
	ptr += sizeof(int);
	
	/* initialize the temporary arrays */
	int *n = new int[l];
	int *m = new int[l];
	int8_t *gh = new int8_t[l];
	double *coeffs = new double[l];
	
	/* load them in */
	for (i=0;i<l;i++) {
		gh[i] = ((int8_t*) ptr)[0];
		ptr += sizeof(int8_t);
	}
	for (i=0;i<l;i++) {
		n[i] = ((int*) ptr)[0];
		ptr += sizeof(int);
	}
	for (i=0;i<l;i++) {
		m[i] = ((int*) ptr)[0];
		ptr += sizeof(int);
	}
	for (i=0;i<l;i++) {
		coeffs[i] = ((double*) ptr)[0];
		ptr += sizeof(double);
	}
	
	/* get n max */
	nmax_ = 0;
	for (i=0;i<l;i++) {
		if (n[i] > nmax_) {
			nmax_ = n[i];
		}
	}
	
	/* calculate the length of the coefficient structure */
	nschc_ = 0;
	for (i=0;i<nmax_;i++) {
		nschc_ += (2 + i);
	}
	
	/* create the structure array */
	schc_ = new struct schmidtcoeffs[nschc_];
	
	/*fill it up */
	p = 0;
	for (i=1;i<=nmax_;i++) {
		for (j=0;j<=i;j++) {
			schc_[p].n = i;
			schc_[p].m = j;
			schc_[p].g = 0.0;
			schc_[p].h = 0.0;
			p++;
		}
	}
	for (i=0;i<l;i++) {
		p = m[i]-1;
		for (j=0;j<n[i];j++) {
			p += (1 + j);
		}
		if (gh[i] == 0) {
			schc_[p].g = coeffs[i];
		} else {
			schc_[p].h = coeffs[i];
		}
	}
			
	/* free the original arrays */
	delete[] n;
	delete[] m;
	delete[] gh;
	delete[] coeffs;
	
}


void Internal::_Schmidt() {
	
	/* create nmax_+1*n+1 array full of coefficients */
	int n, m;
	Snm_ = new double*[nmax_+1];
	for (n=0;n<=nmax_;n++) {
		Snm_[n] = new double[n+1];
	}
	
	/* calculate a bunch of factorials */
	int nfact = 2*nmax_ + 1;
	double facts[nfact];
	facts[0] = 1.0;
	for (n=1;n<nfact;n++) {
		facts[n] = n*facts[n-1];
	}
	
	/* fill the arrays */
	double delta;
	for (n=0;n<=nmax_;n++) {
		for (m=0;m<=n;m++) {
			if (m == 0) {
				delta = 1.0;
			} else { 
				delta = 2.0;
			}
			Snm_[n][m] = sqrt(delta*((facts[n-m]/facts[n+m])));
		}
	}
}

void Internal::_CoeffGrids() {
	
	/* create the grids for g and h */
	int n, m;
	g_ = new double*[nmax_+1];
	h_ = new double*[nmax_+1];
	for (n=0;n<=nmax_;n++) {
		g_[n] = new double[n+1];
		h_[n] = new double[n+1];
		
		/* fill with zeros just in case */
		for (m=0;m<=n;m++) {
			g_[n][m] = 0.0;
			h_[n][m] = 0.0;
			
		}
	}	
	
	/* now we need to fill the arrays in, multiplying by Snm_ */
	for (n=0;n<nschc_;n++) {
		g_[schc_[n].n][schc_[n].m] = schc_[n].g*Snm_[schc_[n].n][schc_[n].m];
		h_[schc_[n].n][schc_[n].m] = schc_[n].h*Snm_[schc_[n].n][schc_[n].m];
	}		
	
		
	
}

void Internal::_Legendre(int l, double *costheta, double *sintheta, 
						int nmax, double ***Pnm, double ***dPnm) {
	
	/* set up the intial few terms */
	int n, m, i;
	for (i=0;i<l;i++) {
		Pnm[0][0][i] = 1.0;
		Pnm[1][0][i] = costheta[i];
		Pnm[1][1][i] = sintheta[i];
		dPnm[0][0][i] = 0.0;
		dPnm[1][0][i] = -sintheta[i];
		dPnm[1][1][i] = costheta[i];
	}
	
	/* now recurse through the rest of them */
	double n21,onenm,nm1;
	for (n=2;n<=nmax;n++) {
		n21 = 2.0*n - 1.0;
		for (m=0;m<=n;m++) {
			if (m < n-1) {
				/* this case is the more complicated one, where we need
				 * two previous polynomials to calculate the next */
				onenm = 1.0/(n-m);
				nm1 = (n + m - 1.0);
				for (i=0;i<l;i++) {
					Pnm[n][m][i] = onenm*(n21*costheta[i]*Pnm[n-1][m][i] - nm1*Pnm[n-2][m][i]);
					dPnm[n][m][i] = onenm*(n21*(costheta[i]*dPnm[n-1][m][i] - sintheta[i]*Pnm[n-1][m][i]) - nm1*dPnm[n-2][m][i]);
				}				
			} else { 
				/* this case only requires one previous polynomial */
				for (i=0;i<l;i++) {
					Pnm[n][m][i] = n21*sintheta[i]*Pnm[n-1][m-1][i];
					dPnm[n][m][i] = n21*(costheta[i]*Pnm[n-1][m-1][i] + sintheta[i]*dPnm[n-1][m-1][i]);
				}
				
				
			}
		}
	}
			
	
}

/* try making a scalar version of this to remove new/delete allocation*/
void Internal::_SphHarm(int l, double *r, double *t, double *p,
					int MaxDeg, double *Br, double *Bt, double *Bp) {
	
	/* set the maximum degree of the model to use */
	int nmax = nmax_;
	if ((MaxDeg > 0) && (MaxDeg < nmax_)) {
		nmax = MaxDeg;
	}
	
	/* create arrays for the Legendre polynomials */
	int n, m, i;
	double ***Pnm = new double**[nmax+1];
	double ***dPnm = new double**[nmax+1];
	for (n=0;n<=nmax;n++) {
		Pnm[n] = new double*[n+1];
		dPnm[n] = new double*[n+1];
		
		for (m=0;m<=n;m++) {
			Pnm[n][m] = new double[l];
			dPnm[n][m] = new double[l];
		}
	}	
	
	/* create some arrays to be used in the field calculation */
	double *r1 = new double[l];
	double *C = new double[l];
	double *cost = new double[l];
	double *sint = new double[l];
	double *sint1 = new double[l];
	for (i=0;i<l;i++) {
		r1[i] = 1.0/r[i];
		C[i] = r1[i]*r1[i];
		cost[i] = cos(t[i]);
		sint[i] = sin(t[i]);
		if (sint[i] == 0.0) {
			sint1[i] = 0.0;
		} else {
			sint1[i] = 1.0/sint[i];
		}
	}
	double **cosmp = new double*[nmax+1];
	double **sinmp = new double*[nmax+1];
	for (m=0;m<=nmax;m++) {
		cosmp[m] = new double[l];
		sinmp[m] = new double[l];
		if (m == 0) {
			for (i=0;i<l;i++) {
				cosmp[0][i] = 1.0;
				sinmp[0][i] = 0.0;
			}
		} else {
			for (i=0;i<l;i++) {
				cosmp[m][i] = cos(((double) m)*p[i]);
				sinmp[m][i] = sin(((double) m)*p[i]);
			}
		}
	}
	double *sumr = new double[l];
	double *sumt = new double[l];
	double *sump = new double[l];
	
	
	/* calculate the Legendre polynomials */
	_Legendre(l,cost,sint,nmax,Pnm,dPnm);
	
	/* set B components to 0 */
	for (i=0;i<l;i++) {
		Br[i] = 0.0;
		Bt[i] = 0.0;
		Bp[i] = 0.0;
	}
	
	/* now start summing stuff up */
	for (n=1;n<=nmax;n++) {
		/* zero the sum arrays and update the C parameter */
		for (i=0;i<l;i++) {
			C[i] = C[i]*r1[i];
			sumr[i] = 0.0;
			sumt[i] = 0.0;
			sump[i] = 0.0;
		}
		
		/* start summing stuff */
		for (m=0;m<=n;m++) {
			for (i=0;i<l;i++) {
				sumr[i] += Pnm[n][m][i]*(g_[n][m]*cosmp[m][i] + h_[n][m]*sinmp[m][i]);
				sumt[i] += dPnm[n][m][i]*(g_[n][m]*cosmp[m][i] + h_[n][m]*sinmp[m][i]);
				sump[i] += ((double) m)*Pnm[n][m][i]*(h_[n][m]*cosmp[m][i] - g_[n][m]*sinmp[m][i]);
			}
		}
		
		/* now calculate B */
		for (i=0;i<l;i++) {
			Br[i] += C[i]*(n+1)*sumr[i];
			Bt[i] += -C[i]*sumt[i];
			Bp[i] += -C[i]*sump[i];
		}
		
	}
	
	/* finally multiply by 1/sintheta */
	for (i=0;i<l;i++) {
		Bp[i] = sint1[i]*Bp[i];
	}
	
	
	/* delete the arrays */
	for (n=0;n<=nmax;n++) {
		for (m=0;m<=n;m++) {
			delete[] Pnm[n][m];
			delete[] dPnm[n][m];
		}
		delete[] Pnm[n];
		delete[] dPnm[n];
		
		delete[] cosmp[n];
		delete[] sinmp[n];
	}		
	delete[] Pnm;
	delete[] dPnm;
	
	delete[] cosmp;
	delete[] sinmp;
	
	delete[] r1;
	delete[] C;
	delete[] cost;
	delete[] sint;
	delete[] sint1;
	
	delete[] sumr;
	delete[] sumt;
	delete[] sump;
					
}

void Internal::_Cart2Pol(int l, double *x, double *y, double *z,
						double *r, double *t, double *p) {
	
	int i;
	double pi2 = M_PI*2;
	for (i=0;i<l;i++) {
		r[i] = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
		t[i] = acos(z[i]/r[i]);
		p[i] = fmod(atan2(y[i],x[i]) + pi2,pi2);
	}
}

void Internal::_BPol2BCart(int l, double *t, double *p,
							double *Br, double *Bt, double *Bp,
							double *Bx, double *By, double *Bz) {
	
	int i;
	double cost, cosp, sint ,sinp;
	for (i=0;i<l;i++) {
		cost = cos(t[i]);
		cosp = cos(p[i]);
		sint = sin(t[i]);
		sinp = sin(p[i]);
		Bx[i] = Br[i]*sint*cosp + Bt[i]*cost*cosp - Bp[i]*sinp;
		By[i] = Br[i]*sint*sinp + Bt[i]*cost*sinp + Bp[i]*cosp;
		Bz[i] = Br[i]*cost - Bt[i]*sint;
	}
								
}

void Internal::Field(int l, double *p0, double *p1, double *p2,
					double *B0, double *B1, double *B2) {
	
	/* some IO pointers */
	int i;
	double *r, *t, *p, *Br, *Bt, *Bp;
	
	/* set the input pointers */
	if (!CartIn_) {
		r = p0;
		t = p1;
		p = p2;
	} else { 
		r = new double[l];
		t = new double[l];
		p = new double[l];
		_Cart2Pol(l,p0,p1,p2,r,t,p);
	}
	
	/* set up the output pointers */
	if (!CartOut_) {
		Br = B0;
		Bt = B1;
		Bp = B2;
	} else { 
		Br = new double[l];
		Bt = new double[l];
		Bp = new double[l];		
	}
	
	/* call the model */
	_SphHarm(l,r,t,p,nmax_,Br,Bt,Bp);
	
	/* rotate field vector if needed and delete output arrays */
	if (CartOut_) {
		_BPol2BCart(l,t,p,Br,Bt,Bp,B0,B1,B2);
		delete[] Br;
		delete[] Bt;
		delete[] Bp;
	}
	
	/* delete input arrays */
	if (CartIn_) {
		delete[] r;
		delete[] t;
		delete[] p;
	}
}

void Internal::Field(int l, double *p0, double *p1, double *p2,
					int MaxDeg, double *B0, double *B1, double *B2) {
	
	/* some IO pointers */
	int i;
	double *r, *t, *p, *Br, *Bt, *Bp;
	
	/* set the input pointers */
	if (!CartIn_) {
		r = p0;
		t = p1;
		p = p2;
	} else { 
		r = new double[l];
		t = new double[l];
		p = new double[l];
		_Cart2Pol(l,p0,p1,p2,r,t,p);
	}
	
	/* set up the output pointers */
	if (!CartOut_) {
		Br = B0;
		Bt = B1;
		Bp = B2;
	} else { 
		Br = new double[l];
		Bt = new double[l];
		Bp = new double[l];		
	}
	
	/* call the model */
	_SphHarm(l,r,t,p,MaxDeg,Br,Bt,Bp);
	
	/* rotate field vector if needed and delete output arrays */
	if (CartOut_) {
		_BPol2BCart(l,t,p,Br,Bt,Bp,B0,B1,B2);
		delete[] Br;
		delete[] Bt;
		delete[] Bp;
	}
	
	/* delete input arrays */
	if (CartIn_) {
		delete[] r;
		delete[] t;
		delete[] p;
	}
}

void Internal::Field(	double p0, double p1, double p2,
						double *B0, double *B1, double *B2) {
	
	/* temporary variables*/
	double r, t, p, Br, Bt, Bp;
	
	/* convert input coords (or not) */
	if (!CartIn_) {
		r = p0;
		t = p1;
		p = p2;
	} else { 
		_Cart2Pol(1,&p0,&p1,&p2,&r,&t,&p);
	}
	
	/* call the model */
	_SphHarm(1,&r,&t,&p,nmax_,&Br,&Bt,&Bp);
	
	/* rotate field vector if needed and delete output arrays */
	if (CartOut_) {
		_BPol2BCart(1,&t,&p,&Br,&Bt,&Bp,B0,B1,B2);
	} else {
		B0[0] = Br;
		B1[0] = Bt;
		B2[0] = Bp;
	}
}

void Internal::Field(	double p0, double p1, double p2, int MaxDeg,
						double *B0, double *B1, double *B2) {
	
	/* temporary variables*/
	double r, t, p, Br, Bt, Bp;
	
	/* convert input coords (or not) */
	if (!CartIn_) {
		r = p0;
		t = p1;
		p = p2;
	} else { 
		_Cart2Pol(1,&p0,&p1,&p2,&r,&t,&p);
	}
	
	/* call the model */
	_SphHarm(1,&r,&t,&p,MaxDeg,&Br,&Bt,&Bp);
	
	/* rotate field vector if needed and delete output arrays */
	if (CartOut_) {
		_BPol2BCart(1,&t,&p,&Br,&Bt,&Bp,B0,B1,B2);
	} else {
		B0[0] = Br;
		B1[0] = Bt;
		B2[0] = Bp;
	}
}
