#include "internal.h"

Internal::Internal(const char *model) {
	
	
	/* use the model string to determine which model to load */
	unsigned char *ptr;
	if (strcmp(model,"VIP4") == 0) {
		/* load the VIP 4 model */
		ptr = &_binary_vip4coeffs_bin_start;
	} else if (strcmp(model,"JRM09") == 0) {
		/* load the JRM09 model */
		ptr = &_binary_jrm09coeffs_bin_start;
	} else {
		/* default to VIP4 */
		ptr = &_binary_vip4coeffs_bin_start;
	}
	
	/* read the coeffs into the object */
	_LoadSchmidt(ptr);
	
	
}

Internal::~Internal() {
	
}

void Internal::_LoadSchmidt(unsigned char *ptr){
	
	/* this is the length of each array */
	int l, i, j, p;
	
	/* read the length */
	l = ((int*) ptr)[0];
	ptr += sizeof(int);
	
	/* initialize the temporary arrays */
	int *n = new int[l];
	int *m = new int[l];
	int8_t *gh = new int8_t[l];
	double *coeffs = new double[l];
	
	/* load them in */
	for (i=0;i<l;i++) {
		gh[i] = ((int8_t*) ptr)[0];
		ptr += sizeof(int8_t);
	}
	for (i=0;i<l;i++) {
		n[i] = ((int*) ptr)[0];
		ptr += sizeof(int);
	}
	for (i=0;i<l;i++) {
		m[i] = ((int*) ptr)[0];
		ptr += sizeof(int);
	}
	for (i=0;i<l;i++) {
		coeffs[i] = ((double*) ptr)[0];
		ptr += sizeof(double);
	}
	
	/* get n max */
	nmax_ = 0;
	for (i=0;i<l;i++) {
		if (n[i] > nmax_) {
			nmax_ = n[i];
		}
	}
	
	/* calculate the length of the coefficient structure */
	nschc_ = 0;
	for (i=0;i<nmax_;i+) {
		nschc_ += (2 + i);
	}
	
	/* create the structure array */
	schc_ = new struct schmidtcoeffs[nschc_];
	
	/*fill it up */
	p = 0;
	for (i=1;i<=nmax_;i++) {
		for (j=0;j<=i;j++) {
			schc_[p].n = i;
			schc_[p].m = j;
			schc_[p].g = 0.0;
			schc_[p].h = 0.0;
			p++;
		}
	}
	for (i=0;i<l;i++) {
		p = m[i]-1
		for (j=0;j<n[i];j++) {
			p += (1 + j);
		}
		if (gh[i] == 0) {
			schc_[p].g = coeffs[i];
		} else {
			schc_[p].h = coeffs[i];
		}
	}
			
	/* free the original arrays */
	delete[] n;
	delete[] m;
	delete[] gh;
	delete[] coeffs;
	
}
