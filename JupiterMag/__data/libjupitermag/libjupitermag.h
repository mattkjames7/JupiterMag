#ifndef __LIBJUPITERMAG_H__
#define __LIBJUPITERMAG_H__
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "trace.h"
#include "string.h"
#include "internal/libinternal.h"
#include "con2020/libcon2020.h"


using namespace std;

extern "C" {
	bool TraceField(int n, double *x0, double *y0, double *z0,
					const char *IntFunc, const char *ExtFunc,
					int MaxLen, double MaxStep, double InitStep,
					double MinStep, double ErrMax, double Delta,
					bool Verbose, int TraceDir,
					int *nstep,
					double **x, double **y, double **z,
					double **Bx, double **By, double **Bz,
					double **R, double **S, double **Rnorm, double **FP,
					int nalpha, double *alpha, double *halpha);
}
#endif
