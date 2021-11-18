#ifndef __TRACECLOSESTPOS_H__
#define __TRACECLOSESTPOS_H__
#include <stdio.h>
#include <stdlib.h>
#include "spline/spline.h"
#include <math.h>
#endif

void traceClosestPos(	int n, double *x, double *y, double *z,
						double *bx, double *by, double *bz,
						int n0, double *x0, double *y0, double *z0,
						int n1, double *x1, double *y1, double *z1,
						double *xc0, double *yc0, double *zc0,
						double *xc1, double *yc1, double *zc1 );

void _ClosestPos(	double px, double py, double pz,
					double bx, double by, double bz,
					int n, double *x, double *y, double *z,
					double *xc, double *yc, double *zc );
					
void _Closest4Pos(	double px, double py, double pz,
					double *x, double *y, double *z, int n,
					int *nc, double *cx, double *cy, double *cz);

void _ClosestPosSpline(int nc, double *cx, double *cy, double *cz,
						double *costheta,
						double *xc, double *yc, double *zc);
						
void _ClosestPosLinear(int nc, double *cx, double *cy, double *cz,
						double *costheta,
						double *xc, double *yc, double *zc);
