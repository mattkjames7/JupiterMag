#ifndef __MODEL_H__
#define __MODEL_H__
#include <stdio.h>
#include <stdlib.h>
#include "libinternalfield/libinternalfield/libinternal.h"
#include "con2020/libcon2020.h"
#endif

extern "C" {
	void ModelField(double p0, double p1, double p2, 
					const char *internal, const char *external, 
					bool CartIn, bool CartOut,
					double *B0, double *B1, double *B2);

	void ModelFieldArray(	int n, double *p0, double *p1, double *p2, 
							const char *internal, const char *external, 
							bool CartIn, bool CartOut,
							double *B0, double *B1, double *B2);
}
