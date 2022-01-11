#include "libjupitermag.h"


bool TraceField(int n, double *x0, double *y0, double *z0,
				const char *IntFunc, const char *ExtFunc,
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
	//if (strcmp(IntFunc,"JRM09") == 0) {
		//Funcs.push_back(JRM09Field);
	//} else if (strcmp(IntFunc,"VIP4") == 0) {
		//Funcs.push_back(VIP4Field);
	//} else if (strcmp(IntFunc,"GSFC13EV") == 0) {
		//Funcs.push_back(GSFC13EVField);
	//} else if (strcmp(IntFunc,"GSFC15EV") == 0) {
		//Funcs.push_back(GSFC15EVField);
	//} else if (strcmp(IntFunc,"GSFC15EVS") == 0) {
		//Funcs.push_back(GSFC15EVSField);
	//} else if (strcmp(IntFunc,"ISAAC") == 0) {
		//Funcs.push_back(ISAACField);
	//} else if (strcmp(IntFunc,"JPL15EV") == 0) {
		//Funcs.push_back(JPL15EVField);
	//} else if (strcmp(IntFunc,"JPL15EVS") == 0) {
		//Funcs.push_back(JPL15EVSField);
	//} else if (strcmp(IntFunc,"O4") == 0) {
		//Funcs.push_back(O4Field);
	//} else if (strcmp(IntFunc,"O6") == 0) {
		//Funcs.push_back(O6Field);
	//} else if (strcmp(IntFunc,"P11A") == 0) {
		//Funcs.push_back(P11AField);
	//} else if (strcmp(IntFunc,"SHA") == 0) {
		//Funcs.push_back(SHAField);
	//} else if (strcmp(IntFunc,"U17EV") == 0) {
		//Funcs.push_back(U17EVField);
	//} else if (strcmp(IntFunc,"V117EV") == 0) {
		//Funcs.push_back(V117EVField);
	//} else if (strcmp(IntFunc,"VIPAL") == 0) {
		//Funcs.push_back(VIPALField);
	//} else if (strcmp(IntFunc,"VIT4") == 0) {
		//Funcs.push_back(VIT4Field);
	//}
	Funcs.push_back(getModelFieldPtr(IntFunc));
	

	/* external model */
	if (strcmp(ExtFunc,"Con2020") == 0) {
		Funcs.push_back(Con2020Field);
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
