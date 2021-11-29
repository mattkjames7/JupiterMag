#include "traceclosestpos.h"


void traceClosestPos(	int n, double *x, double *y, double *z,
						double *bx, double *by, double *bz,
						int n0, double *x0, double *y0, double *z0,
						int n1, double *x1, double *y1, double *z1,
						double *xc0, double *yc0, double *zc0,
						double *xc1, double *yc1, double *zc1 ) {
							
	/* return an array of posisions along two nearby field lines which
	 * are the closest points to each element of the original field line */
	
	
	/* calculate the closest position for each step */
	int i;
	double d0, d1;
	for (i=0;i<n;i++) {
		_ClosestPos(	x[i],y[i],z[i],
						bx[i],by[i],bz[i],
						n0,x0,y0,z0,
						&xc0[i],&yc0[i],&zc0[i]);
		_ClosestPos(	x[i],y[i],z[i],
						bx[i],by[i],bz[i],
						n1,x1,y1,z1,
						&xc1[i],&yc1[i],&zc1[i]);
		d0 = sqrt(pow(x[i]-xc0[i],2) + pow(y[i]-yc0[i],2) + pow(z[i]-zc0[i],2));
		d1 = sqrt(pow(x[i]-xc1[i],2) + pow(y[i]-yc1[i],2) + pow(z[i]-zc1[i],2));
		printf("%d %f %f\n",i,d0,d1);
	}
							
}

void _ClosestPos(	double px, double py, double pz,
					double bx, double by, double bz,
					int n, double *x, double *y, double *z,
					double *xc, double *yc, double *zc ) {

	
	/* find the four(hopefully) closest values */
	int nc;
	double cx[4], cy[4], cz[4];
	_Closest4Pos(px,py,pz,n,x,y,z,&nc,cx,cy,cz);
	
	/* get the unit vector for B*/
	double B = sqrt(bx*bx + by*by + bz*bz);
	double uBx = bx/B;
	double uBy = by/B;
	double uBz = bz/B;
	
	/* get unit vectors for the vector fromt he current position to
	 * each of the four closest positions on the adjacent field line */
	double udx[4], udy[4], udz[4], D;
	int i;
	for (i=0;i<nc;i++) {
		udx[i] = cx[i] - px;
		udy[i] = cy[i] - py;
		udz[i] = cz[i] - pz;
		D = sqrt(udx[i]*udx[i] + udy[i]*udy[i] + udz[i]*udz[i]);
		udx[i] = udx[i]/D;
		udy[i] = udy[i]/D;
		udz[i] = udz[i]/D;
	}
	
	/* calculate the angle between the unit vectors using the dot product*/
	double costheta[4];
	for (i=0;i<nc;i++) {
		costheta[i] = udx[i]*uBx + udy[i]*uBy + udz[i]*uBz;
		printf("%d %f\n",i,acos(costheta[i])*180/M_PI);
	}
	
	/* interpolate to find approximately where cos(theta) = 0 */
	if (nc == 4) {
		_ClosestPosSpline(nc,cx,cy,cz,costheta,xc,yc,zc);
	} else {
		_ClosestPosLinear(nc,cx,cy,cz,costheta,xc,yc,zc);
	}
}



void _Closest4Pos(	double px, double py, double pz,
					int n, double *x, double *y, double *z,
					int *nc, double *cx, double *cy, double *cz) {
	int i, k=0;
	if (n == 0) {
		nc[0] = 0;
		return;
	}
	
	if (n <= 4) {
		nc[0] = n;
		for (i=0;i<nc[0];i++) {
			cx[i] = x[i];
			cy[i] = y[i];
			cz[i] = z[i];
		}
		return;
	}

	double dx, dy, dz;
	double d;
	double dmin = INFINITY;
	int Imind = -1;

	/*find the closest point */
	for (i=0;i<n;i++) {
		dx = px - x[i];
		dy = py - y[i];
		dz = pz - z[i];
		d = dx*dx + dy*dy + dz*dz;
		if (d < dmin) {
			dmin = d;
			Imind = i;
		}
	}
	printf("I: %d D: %f\n",Imind,sqrt(dmin));

	int I4[4], imn, imx;

	if (Imind == (n-1)) {
		Imind--;
	}

	if (z[Imind] > z[Imind+1]) {
		I4[0] = Imind - 1;
		I4[1] = Imind;
		I4[2] = Imind + 1;
		I4[3] = Imind + 2;
		imn = 0;
		imx = 3;
	} else {
		I4[3] = Imind - 1;
		I4[2] = Imind;
		I4[1] = Imind + 1;
		I4[0] = Imind + 2;
		imn = 3;
		imx = 0;
	}		

	while (I4[imn] < 0) {
		for (i=0;i<4;i++) {
			I4[i]++;
		}
	}

	while (I4[imx] >= n) {
		for (i=0;i<4;i++) {
			I4[i]--;
		}
	}
	nc[0] = 4;

	for (i=0;i<nc[0];i++) {
		cx[i] = x[I4[i]];
		cy[i] = y[I4[i]];
		cz[i] = z[I4[i]];
		//printf("%f ",sqrt(pow(px - cx[i],2) + pow(py - cy[i],2) + pow(pz - cz[i],2)));
	}
	//printf("\n");
	nc[0] = 4;
	
}


void _ClosestPosSpline(int nc, double *cx, double *cy, double *cz,
						double *costheta,
						double *xc, double *yc, double *zc) {
	
	/* create some splines*/
	Spline Sx(nc,costheta,cx);
	Spline Sy(nc,costheta,cy);
	Spline Sz(nc,costheta,cz);
	
	/* interpolate */
	double ct0 = 0.0;
	Sx.Interpolate(1,&ct0,xc);
	Sy.Interpolate(1,&ct0,yc);
	Sz.Interpolate(1,&ct0,zc);
}

void _ClosestPosLinear(int nc, double *cx, double *cy, double *cz,
						double *costheta,
						double *xc, double *yc, double *zc) {
	
	/* pick the closest pair of positions */
	int i, i0=-1, i1=-1;
	for (i=0;i<nc-1;i++) {
		if ((0.0 >= costheta[i]) && (0.0 < costheta[i+1])) {
			i0 = i;
			i1 = i + 1;
		}
	}
	
	if (i0 == -1) {
		if (0.0 < costheta[0]) {
			i0 = 0;
			i1 = 1;
		}
		if (0.0 >= costheta[nc-1]) {
			i0 = nc - 2;
			i1 = nc - 1;
		}
	}
	
	/* interpolate */
	double dct = 0.0 - costheta[i0];
	double m;
	m = (cx[i1] - cx[i0])/(costheta[i1] - costheta[i1]);
	xc[0] = m*dct + cx[i0];
	m = (cy[i1] - cy[i0])/(costheta[i1] - costheta[i1]);
	yc[0] = m*dct + cy[i0];
	m = (cz[i1] - cz[i0])/(costheta[i1] - costheta[i1]);
	zc[0] = m*dct + cz[i0];
	

}
