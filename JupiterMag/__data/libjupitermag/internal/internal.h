#ifndef __INTERNAL_H__
#define __INTERNAL_H__
#include <stdio.h>
#include <stdlib.h>

#endif
using namespace std;

/* here are the pointers to the memory where the coefficients are stored*/
extern unsigned char _binary_vip4coeffs_bin_start;
extern unsigned char _binary_jrm09coeffs_bin_start;

/* This structure will store the Schmidt coefficients */
struct schmidtcoeffs {
	int n;
	int m;
	double g;
	double h;
};

class Internal {
	public:
		Internal(const char*);
		~Internal();
	
		/*these two functions will calculate the field in Cartesian RH 
		 * system III coordinates for a scalar or an array, 
		 * respectively.*/
		void Fields(double,double,double,int,double*,double*,double*);
		void Fielda(int,double*,double*,double*,int,double*,double*,double*);
		
		
	private:
		/*Schmidt coefficients */
		struct schmidtcoeffs *schc_;
		int nschc_;
		double **Snm_;
		int nmax_;
		
		/* these ones will have Snm_ already multiplied */
		double **g_;
		double **h_;
		
		/* functions for initializing the object */
		void _LoadSchmidt(const char*);
		void _Schmidt();
		void _CoeffGrids();

		
	
};
