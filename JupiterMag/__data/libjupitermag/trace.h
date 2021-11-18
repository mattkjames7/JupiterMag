#ifndef __TRACE_H__
#define __TRACE_H__
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include "traceclosestpos.h"
#endif
using namespace std;
/* this will be used for all of the model wrapper functions (configure model first) */
typedef void (*FieldFuncPtr)(double,double,double,double*,double*,double*);



class Trace {
	
	public:
		Trace(vector<FieldFuncPtr>);
		~Trace();
		
		void InputPos(int,double*,double*,double*);
		void SetTraceCFG(int, double,double,double,double,double,bool,int);
		void SetTraceCFG();
		
		void SetAlpha(int,double*,double);
	
		/* tracing */
		void TraceField(int*,double**,double**,double**,double**,double**,double**,double**);
		void TraceField();
		void StepVector(double,double,double,double,double*,double*,double*);
		bool ContinueTrace(double,double,double,double*);
		void Step(double,double,double,double*,double*,double*,double*,double*,double*,double*);
		void ReverseElements(int, double*);
		void RKMTrace(	double,double,double,int*,double*,
						double*,double*,double*,double*,double*,double*)

		/* get a single field vector */
		void Field(double,double,double,double*,double*,double*);

		/* calculate trace distance,R,Rnorm */
		void CalculateTraceDist(double**);
		void CalculateTraceDist();
		void _CalculateTraceDist();
		void CalculateTraceRnorm(double**);
		void CalculateTraceRnorm();
		void _CalculateTraceRnorm();
	
		/* Calculate footprints */
		void CalculateTraceFP(double**);
		void CalculateTraceFP();
		void _CalculateTraceFP();
		
		/* calculate halpha */
		void CalculateHalpha();
		void CalculateHalpha(double*);
		void CalculateHalpha(double***);
		void CalculateHalpha(double*,double***);
	
		/* return things*/
		void GetTraceNstep(int*);
		void GetTrace(double**,double**,double**);
		void GetTrace(double**,double**,double**,double**,double**,double**);
		void GetTraceDist(double**);
		void GetTraceR(double**);
		void GetTraceRnorm(double**);
		void GetTraceFootprints(double**);
		void GetTraceHalpha(double*);	/* python will use this */
		void GetTraceHalpha(double***); /* no idea how to link this to python*/
		
		Trace TracePosition(int,double,double,double);

		/* input coords */
		int n_;
		double *x0_, *y0_, *z0_;  
		int *Date_;
		float *ut_;

		/* trace params */
		int MaxLen_;
		double MaxStep_, MinStep_, InitStep_;
		bool Verbose_;
		int TraceDir_;
		double ErrMax_;
		
		/* trace coords */
		int *nstep_;
		double **x_, **y_, **z_;
	
		/* trace fields */
		double **bx_, **by_, **bz_;

		/* trace end points */
		double *xfn_, *yfn_, *zfn_;
		double *xfs_, *yfs_, *zfs_;
		double *xfe_, *yfe_, *zfe_;
	
	private:
		/* this is the number of field contributions */
		int nf_;
		vector<FieldFuncPtr> Funcs_;

		/* booleans to tell the object what has been done */
		bool inputPos_;
		bool tracedField_,allocTrace_;
		bool hasFootprints_,allocFootprints_;
		bool allocEqFP_;
		bool hasDist_,allocDist_;
		bool hasRnorm_,allocRnorm_;
		bool hasHalpha_,allocHalpha_, allocHalpha3D_;
		bool allocAlpha_;

		

	
		/* field length, R, Rnorm, Halpha, Footprints */
		int nalpha_;
		double *alpha0_, *alpha1_;
		double Delta_;
		double **S_;
		double **R_;
		double **Rnorm_;
		double *Halpha_;
		double ***Halpha3D_;
		double **FP_;
		
		/* model */
		const char *Model_;
		ModelFuncPtr ModelFunc_;
	
		/* hidden trace functions */
		void _TraceField();

		/* halpha functions */
		bool _CheckHalpha();
		void _CalculateHalpha();
		void _CalculateTraceHalpha(int,int,double*);
		void _CalculateHalphaStartPoints(int i, int j,
							double *xe0, double *ye0, double *ze0,
							double *xe1, double *ye1, double *ze1);
	
};
