#ifndef __INTERNAL_H__
#define __INTERNAL_H__
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#endif
using namespace std;

/* here are the pointers to the memory where the coefficients are stored*/
extern unsigned char _binary_vip4coeffs_bin_start;
extern unsigned char _binary_jrm09coeffs_bin_start;
extern unsigned char _binary_gsfc13evcoeffs_bin_start;
extern unsigned char _binary_gsfc15evcoeffs_bin_start;
extern unsigned char _binary_gsfc15evscoeffs_bin_start;
extern unsigned char _binary_isaaccoeffs_bin_start;
extern unsigned char _binary_jpl15evcoeffs_bin_start;
extern unsigned char _binary_jpl15evscoeffs_bin_start;
extern unsigned char _binary_o4coeffs_bin_start;
extern unsigned char _binary_o6coeffs_bin_start;
extern unsigned char _binary_p11acoeffs_bin_start;
extern unsigned char _binary_shacoeffs_bin_start;
extern unsigned char _binary_u17evcoeffs_bin_start;
extern unsigned char _binary_v117evcoeffs_bin_start;
extern unsigned char _binary_vipalcoeffs_bin_start;
extern unsigned char _binary_vit4coeffs_bin_start;

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
		Internal(unsigned char *);
		~Internal();
	
		/*these two functions will calculate the field in Cartesian RH 
		 * system III coordinates.*/
		void Field(int,double*,double*,double*,double*,double*,double*);
		void Field(double,double,double,double*,double*,double*);
		
		/* set model parameters*/
		void SetCartIn(bool);
		void SetCartOut(bool);
		bool GetCartIn();
		bool GetCartOut();
		
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
		void _LoadSchmidt(unsigned char*);
		void _Schmidt();
		void _CoeffGrids();

		/* This function will calculate the Legendre polynomials */
		void _Legendre(int,double*,double*,double***,double***);
		
		/* this function will calculate the magnetic field components in
		 * spherical polar coordinates */
		void _SphHarm(int,double*,double*,double*,double*,double*,double*);
		
		/* coordinate/field vector rotation */
		bool CartIn_;
		bool CartOut_;
		void _Cart2Pol(int,double*,double*,double*,double*,double*,double*);
		void _BPol2BCart(int,double*,double*,double*,double*,double*,double*,double*,double*);
	
};
