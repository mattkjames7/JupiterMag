#ifndef __BESSEL_H__
#define __BESSEL_H__
#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>

#endif
using namespace std;

void j0(int n, double *x, double *j);
void j1(int n, double *x, double *j);

void j0(int n, double *x, double multx, double *j);
void j1(int n, double *x, double multx, double *j);
