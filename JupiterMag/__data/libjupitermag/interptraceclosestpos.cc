#include "interptraceclosestpos.h"

void interptraceClosestPos(	int n, double *x, double *y, double *z,
						double *bx, double *by, double *bz,
						int n0, double *x0, double *y0, double *z0, double *s0,
						int n1, double *x1, double *y1, double *z1, double *s1
						double *xc0, double *yc0, double *zc0,
						double *xc1, double *yc1, double *zc1 ) {
							
	/* return an array of posisions along two nearby field lines which
	 * are the closest points to each element of the original field line */
	
	/* get a couple of splines (one for each field line) */
	Spline Sx0(n0,s0,x0);
	Spline Sy0(n0,s0,y0);
	Spline Sz0(n0,s0,z0);
	Spline Sx1(n1,s1,x1);
	Spline Sy1(n1,s1,y1);
	Spline Sz1(n1,s1,z1);
	
	

	/* find the closest position along the splines for each position */
	int i;
	double s0_0, s1_0;
	for (i=0;i<n;i++) {
		/* find the closest trace position to start the optimization at */
		s0_0 = ClosestS(x[i],y[i],z[i],n0,x0,y0,z0,s0);
		s1_0 = ClosestS(x[i],y[i],z[i],n1,x1,y1,z1,s1);
		
		/* Trace 0 */
		OptimisePos(x[i],y[i],z[i],bx[i],by[i],bz[i],s0,Sx0,Sy0,Sz0,&xc0[i],&yc0[i],&zc0[i]);
		
		/* Trace 1 */
		OptimisePos(x[i],y[i],z[i],bx[i],by[i],bz[i],s1,Sx1,Sy1,Sz1,&xc1[i],&yc1[i],&zc1[i]);
	}
							
}

double ClosestS(double x, double y, double z,
				int nt, double *xt, double *yt, double *zt,
				double *st) {

	int i, imin;
	double dx, dy, dz, d, dmin = INF;
	for (i=0;i<nt,i++) {
		dx = x - xt[i];
		dy = y - yt[i];
		dz = z - zt[i];
		d = sqrt(dx*dx + dy*dy + dz*dz);
		if (d < dmin) {
			imin = i;
			dmin = d;
		}
	}
	return st[imin];
					
}

void OptimizePos(	double x, double y, double z,
					double bx, double by, double bz,
					double s0, 
					Spline Sx, Spline Sy, Spline Sz,
					double *xc, double *yc, double *zc) {
						
						
}