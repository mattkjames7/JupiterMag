#include "test.h"

int main() {
	
	int n = 1;
	double x0 = 5.0;
	double y0 = 0.0;
	double z0 = 0.0;
	int nalpha = 1;
	double alpha = 0.0;
	
	printf("Create field function vectors\n");
	/* test function for debugging the trace */
	vector<FieldFuncPtr> Funcs;

	/* internal model */
	Funcs.push_back(jrm09Field);

	/* external model */
	Funcs.push_back(Con2020Field);

	/* initialise the trace object */
	printf("Create Trace object\n");
	Trace T(Funcs);

	/* add the starting posiutions fo the traces */
	printf("Add starting position\n");
	T.InputPos(n,&x0,&y0,&z0);

	/* configure the trace parameters */
	printf("Set the trace parameters \n");
	T.SetTraceCFG();

	/* set up the alpha calculation */
	printf("Initialize alpha\n");
	T.SetAlpha(nalpha,&alpha);
	

	/* Trace */
	printf("Trace\n");
	T.TraceField();

	/* trace distance, footprints, Rnorm */
	printf("Footprints etc...\n");
	T.CalculateTraceDist();
	T.CalculateTraceFP();
	T.CalculateTraceRnorm();
	

	/* halpha */
	printf("H_alpha\n");
	T.CalculateHalpha();
	
	printf("Trace s and h_alpha\n");
	int i;
	for (i=0;i<T.nstep_[0];i++) {
		printf("%10.4f %10.4f\n",T.S_[0][i],T.Halpha3D_[0][0][i]);
	}

	
}
