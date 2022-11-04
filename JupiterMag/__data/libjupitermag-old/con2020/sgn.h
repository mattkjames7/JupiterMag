#ifndef __SGN_H__
#define __SGN_H__
#include <stdio.h>
#include <stdlib.h>

#endif
using namespace std;

template <typename T> T sgn(T x) {
	return (x > 0) - (x < 0);
}
