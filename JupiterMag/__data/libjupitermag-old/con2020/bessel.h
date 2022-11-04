#ifndef __BESSEL_H__
#define __BESSEL_H__
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include "polyeval.h"

/***********************************************************************
 * NAME : j0(x)
 * 
 * DESCRIPTION : Fucntion to calculate an estimate of the Bessel 
 * function j0 using code based on the Cephes C library
 * (Cephes Mathematical Functions Library, http://www.netlib.org/cephes/).
 * 
 * INPUTS :
 * 		double 	x	position to calculate J0 at.
 * 
 * RETURNS :
 * 		double j	j0 function evaluated at x.
 * 
 * ********************************************************************/
double j0(double x);

/***********************************************************************
 * NAME : j1(x)
 * 
 * DESCRIPTION : Fucntion to calculate an estimate of the Bessel 
 * function j1 using code based on the Cephes C library
 * (Cephes Mathematical Functions Library, http://www.netlib.org/cephes/).
 * 
 * INPUTS :
 * 		double 	x	position to calculate J1 at.
 * 
 * RETURNS :
 * 		double j	j1 function evaluated at x.
 * 
 * ********************************************************************/
double j1(double x);

/***********************************************************************
 * NAME : j0(n,x,j)
 * 
 * DESCRIPTION : Fucntion to calculate an estimate of the Bessel 
 * function j0 using code based on the Cephes C library
 * (Cephes Mathematical Functions Library, http://www.netlib.org/cephes/).
 * 
 * INPUTS :
 * 		int 	n	Number of elements in x
 * 		double 	*x	position to calculate J0 at.
 * 
 * OUTPUTS :
 * 		double *j	j0 function evaluated at x.
 * 
 * ********************************************************************/
void j0(int n, double *x, double *j);

/***********************************************************************
 * NAME : j1(n,x,j)
 * 
 * DESCRIPTION : Fucntion to calculate an estimate of the Bessel 
 * function j1 using code based on the Cephes C library
 * (Cephes Mathematical Functions Library, http://www.netlib.org/cephes/).
 * 
 * INPUTS :
 * 		int 	n	Number of elements in x
 * 		double 	*x	position to calculate J1 at.
 * 
 * OUTPUTS :
 * 		double *j	j1 function evaluated at x.
 * 
 * ********************************************************************/
void j1(int n, double *x, double *j);

/***********************************************************************
 * NAME : j0(n,x,multx,j)
 * 
 * DESCRIPTION : Fucntion to calculate an estimate of the Bessel 
 * function j0 using code based on the Cephes C library
 * (Cephes Mathematical Functions Library, http://www.netlib.org/cephes/).
 * 
 * INPUTS :
 * 		int 	n	Number of elements in x
 * 		double 	*x	position to calculate J0(x*multx) at.
 * 		double multx	Constant to multiply x by
 * 
 * OUTPUTS :
 * 		double *j	j0 function evaluated at x*multx.
 * 
 * ********************************************************************/
void j0(int n, double *x, double multx, double *j);

/***********************************************************************
 * NAME : j1(n,x,multx,j)
 * 
 * DESCRIPTION : Fucntion to calculate an estimate of the Bessel 
 * function j1 using code based on the Cephes C library
 * (Cephes Mathematical Functions Library, http://www.netlib.org/cephes/).
 * 
 * INPUTS :
 * 		int 	n	Number of elements in x
 * 		double 	*x	position to calculate J1(x*multx) at.
 * 		double multx	Constant to multiply x by
 * 
 * OUTPUTS :
 * 		double *j	j1 function evaluated at x*multx.
 * 
 * ********************************************************************/
void j1(int n, double *x, double multx, double *j);



#endif
