#ifndef __CLIP_H__
#define __CLIP_H__
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#endif
using namespace std;



template <typename T> T clip(T x, T mn, T mx) {
	return min(mx,max(x,mn));
}
