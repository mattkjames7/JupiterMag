#include "test.h"

int main() {
	
	int n = 1;
	double x0 = 5.0;
	double y0 = 0.0;
	double z0 = 0.0;
	int nalpha = 1;
	double alpha = 0.0;
	
	/* test function for debugging the trace */
	vector<FieldFuncPtr> Funcs;

	/* internal model */
	Funcs.push_back(JRM09Field);

	/* external model */
	Funcs.push_back(Con2020Field);

	/* initialise the trace object */
	Trace T(Funcs);

	/* add the starting posiutions fo the traces */
	T.InputPos(n,&x0,&y0,&z0);

	/* configure the trace parameters */
	T.SetTraceCFG();

	/* set up the alpha calculation */
	T.SetAlpha(nalpha,&alpha);
	

	/* Trace */
	T.TraceField();

	/* trace distance, footprints, Rnorm */
	T.CalculateTraceDist();
	T.CalculateTraceFP();
	T.CalculateTraceRnorm();
	

	/* halpha */
	T.CalculateHalpha();
	
	printf("Trace s and h_alpha\n");
	int i;
	for (i=0;i<T.nstep_[0];i++) {
		printf("%10.4f %10.4f\n",T.S_[0][i],T.Halpha3D_[0][0][i]);
	}

	
}
