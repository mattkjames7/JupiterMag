#ifndef __CON2020_H__
#define __CON2020_H__
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "bessel.h"
#include "sgn.h"
#define deg2rad = M_PI/180.0;
#endif
using namespace std;


class Con2020 {
	public:
		/* constructors */
		Con2020();
		Con2020(double,double,double,double,double,double,double,const char*);
	
		/* destructor */
		~Con2020();
		
		/* these functions will be used to set the equations used, if
		 * they need to be changed post-initialisation */
		void UseEdwardsEqs(bool);
		void SetEqType(const char*);
		
		/* This function will be used to call the model */
		Field(int,double*,double*,double*,double*,double*,double*,bool,bool);
		
	private:
		/* model parameters */
		double mui_,irho_,r0_,r1_,d_,xt_,xp_,dipshift_,diptilt_;
		const char *eqtype_;
		bool Edwards_;
		
		/* Bessel function arrays - arrays prefixed with r and z are
		 * to be used for integrals which calcualte Brho and Bz,
		 * respectively */
		int *rnbes_;			/* number of elements for each Bessel function (rho)*/
		int *znbes_;			/* same as above for z */
		double **rlambda_;/* Lambda array to integrate over rho*/
		double **zlambda_;/* Lambda array to integrate over z*/
		double **rj0_lambda_r0_; /* j0(lambda*r0) */
		double **rj1_lambda_rho_;/* j1(lambda*rho) */
		double **zj0_lambda_r0_; /* j0(lambda*r0) */
		double **zj0_lambda_rho_;/* j0(lambda*rho) */

		
		/* arrays to multiply be stuff to be integrated */
		/* these arrays will store the parts of equations 14, 15, 17 
		 * and 18 of Connerny 1981 which only need to be calculated once*/
		double **Eq14_;		/* j0(lambda*r0)*sinh(lamba*d)/lambda */
		double **Eq15_;     /* j0(lambda*r0)*sinh(lamba*d)/lambda */
		double **Eq17_;     /* j0(lambda*r0)*exp(-lamba*d)/lambda */
		double **Eq18_;     /* j0(lambda*r0)/lambda */
		double **ExpLambdaD_;
		

		/* integration step sizes */
		double dlambda_ = 1e-4;
		double dlambda_brho_ = 1e-4;
		double dlambda_bz_ = 5e-5;
		
		/* Arrays containing maximum lambda values */
		double rlmx_array_[6];
		double zlmx_array_[6];


		/* coordinate conversions for positions */
		void _SysIII2Mag(int,double*,double*,double*,double*,double*,double*,double*);
		void _PolSysIII2Mag(int,double*,double*,double*,double*,double*,double*,double*);
		
		
		/* coordinate conversion for magnetic field vector */
		
		/* Azimuthal field */
		_AzimuthalField(int,double*,double*);
		
};
