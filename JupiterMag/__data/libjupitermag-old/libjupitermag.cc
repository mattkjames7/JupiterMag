#include "libjupitermag.h"


bool TraceField(int n, double *x0, double *y0, double *z0,
				const char *IntFunc, int nExt, char **ExtFunc,
				int MaxLen, double MaxStep, double InitStep,
				double MinStep, double ErrMax, double Delta,
				bool Verbose, int TraceDir,
				int *nstep,
				double **x, double **y, double **z,
				double **Bx, double **By, double **Bz,
				double **R, double **S, double **Rnorm, double **FP,
				int nalpha, double *alpha, double *halpha) {
	
	/* before calling this wrapper function, any field models used should
	 * be configured. This function will not do any of that, so strange
	 * things could happen. Make sure that all models are Cartesian in 
	 * and out! */
	vector<FieldFuncPtr> Funcs;

	/* internal model */
	Funcs.push_back(getModelFieldPtr(IntFunc));

	/* external model */
	int i;
	for (i=0;i<nExt;i++) {
		if (strcmp(ExtFunc[i],"Con2020") == 0) {
			Funcs.push_back(Con2020Field);
		}
	}

	/* if there are no functions then return */
	if (Funcs.size() == 0) {
		printf("No valid model functions provided\n");
		return false;
	}

	/* initialise the trace object */
	Trace T(Funcs);

	/* add the starting posiutions fo the traces */
	T.InputPos(n,x0,y0,z0);

	/* configure the trace parameters */
	T.SetTraceCFG(MaxLen,MaxStep,InitStep,MinStep,ErrMax,Delta,Verbose,TraceDir);

	/* set up the alpha calculation */
	if (nalpha > 0) {
		T.SetAlpha(nalpha,alpha);
	}

	/* Trace */
	T.TraceField(nstep,x,y,z,R,Bx,By,Bz);

	/* trace distance, footprints, Rnorm */
	if (TraceDir == 0) {
		T.CalculateTraceDist(S);
		T.CalculateTraceFP(FP);
		T.CalculateTraceRnorm(Rnorm);
	}

	/* halpha */
	if ((nalpha > 0) && (TraceDir == 0)) {
		T.CalculateHalpha(halpha);
	}

	return true;
}
