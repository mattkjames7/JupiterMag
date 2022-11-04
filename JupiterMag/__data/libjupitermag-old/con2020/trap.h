#ifndef __TRAP_H__
#define __TRAP_H__
#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#endif


double trap(int n, double *x, double *y);
double trapc(int n, double dx, double *y);
