#ifndef __CON2020_H__
#define __CON2020_H__
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "bessel.h"
#include "sgn.h"
#include "clip.h"
#include "trap.h"
#define deg2rad M_PI/180.0;
#endif
using namespace std;

/* function pointer for input conversion */
class Con2020; /*this is needed for the pointer below */ 
typedef void (Con2020::*InputConvFunc)(int,double*,double*,double*,
						double*,double*,double*,double*,double*,
						double*,double*,double*,double*);
/* Output conversion */
typedef void (Con2020::*OutputConvFunc)(int,double*,double*,double*,
				double*,double*,double*,double*,
				double*,double*,double*,
				double*,double*,double*);

/* Model function */
typedef void (Con2020::*ModelFunc)(double,double,double,double*,double*,double*);

/* analytical approximation equations */
typedef void (Con2020::*Approx)(double,double,double,double,double,double*,double*);

class Con2020 {
	public:
		/* constructors */
		Con2020();
		Con2020(double,double,double,double,double,double,double,const char*,bool,bool,bool,bool);
	
		/* destructor */
		~Con2020();
		
		/* these functions will be used to set the equations used, if
		 * they need to be changed post-initialisation */
		void SetEdwardsEqs(bool);
		void SetEqType(const char*);
		void SetCurrentDensity(double);
		void SetRadCurrentDensity(double);
		void SetR0(double);
		void SetR1(double);
		void SetCSHalfThickness(double);
		void SetCSTilt(double);
		void SetCSTiltAzimuth(double);
		void SetErrCheck(bool);
		void SetCartIn(bool);
		void SetCartOut(bool);
		
		/* these mamber functions will be the "getter" version of the
		 * above setters */
		bool GetEdwardsEqs();
		void GetEqType(char*);
		double GetCurrentDensity();
		double GetRadCurrentDensity();
		double GetR0();
		double GetR1();
		double GetCSHalfThickness();
		double GetCSTilt();
		double GetCSTiltAzimuth();
		bool GetErrCheck();
		bool GetCartIn();
		bool GetCartOut();
		
		/* This function will be used to call the model, it is overloaded
		 * so that we have one for arrays, one for scalars */
		void Field(int,double*,double*,double*,double*,double*,double*);
		void Field(double,double,double,double*,double*,double*);
		
	private:
		/* model parameters */
		double mui_,irho_,r0_,r1_,d_,xt_,xp_,disctilt_,discshift_;
		double r0sq_, r1sq_;
		double cosxp_,sinxp_,cosxt_,sinxt_;
		char eqtype_[9];
		bool Edwards_, ErrChk_;
		bool CartIn_,CartOut_;
		
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
		static constexpr double dlambda_ = 1e-4;
		static constexpr double dlambda_brho_ = 1e-4;
		static constexpr double dlambda_bz_ = 5e-5;
		
		/* Arrays containing maximum lambda values */
		double rlmx_array_[6];
		double zlmx_array_[6];


		/* coordinate conversions for positions */
		InputConvFunc _ConvInput;
		void _SysIII2Mag(int,double*,double*,double*,
						double*,double*,double*,double*,double*,
						double*,double*,double*,double*);
		void _PolSysIII2Mag(int,double*,double*,double*,
						double*,double*,double*,double*,double*,
						double*,double*,double*,double*);
		
		
		/* coordinate conversion for magnetic field vector */
		OutputConvFunc _ConvOutput;
		void _BMag2SysIII(int,double*,double*,double*,
							double*,double*,double*,double*,
							double*,double*,double*,
							double*,double*,double*);
		void _BMag2PolSysIII(int,double*,double*,double*,
							double*,double*,double*,double*,
							double*,double*,double*,
							double*,double*,double*);	

		/* Functions to update function pointers */
		void _SetIOFunctions();
		void _SetModelFunctions();
		ModelFunc _Model;
		
							
		/* Azimuthal field */
		void _AzimuthalField(int,double*,double*,double*,double*);
		void _AzimuthalField(double,double,double,double*);
		
		/* analytic equations */
		void _Analytic(double,double,double,double*,double*,double*);
		void _SolveAnalytic(int,double*,double*,double,double*,double*);
		Approx _LargeRho;
		Approx _SmallRho;
		void _AnalyticInner(double,double,double*,double*);
		void _AnalyticOuter(double,double,double*,double*);
		void _LargeRhoConnerney(double,double,double,double,double,double*,double*);
		void _SmallRhoConnerney(double,double,double,double,double,double*,double*);
		void _LargeRhoEdwards(double,double,double,double,double,double*,double*);
		void _SmallRhoEdwards(double,double,double,double,double,double*,double*);
		
		/* integral-related functions */
		void _Integral(double,double,double,double*,double*,double*);
		void _IntegralInner(double, double, double,	double*, double*);
		void _InitIntegrals();
		void _RecalcIntegrals();
		void _DeleteIntegrals();
		void _IntegralChecks(int,double*,int*,int[]);
		void _IntegralCheck(double,int*);
		void _SolveIntegral(int,double*,double*,double*,double*,double*);
		void _IntegrateEq14(int,double,double,double,double*);
		void _IntegrateEq15(int,double,double,double*);
		void _IntegrateEq17(int,double,double,double*);
		void _IntegrateEq18(int,double,double,double*);

		/* hybrid */
		void _Hybrid(double,double,double,double*,double*,double*);
};
