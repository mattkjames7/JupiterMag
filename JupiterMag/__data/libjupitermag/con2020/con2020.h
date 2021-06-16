#ifndef __CON2020_H__
#define __CON2020_H__
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "bessel.h"
#define deg2rad = M_PI/180.0;
#endif
using namespace std;


class Con2020 {
	public:
		/* constructors */
		Con2020();
		Con2020(double,double,double,double,double,double,double,const char*);
	
		/* destructor */
		~Con2020();
		
		/* these functions will be used to set the equations used, if
		 * they need to be changed post-initialisation */
		void UseEdwardsEqs(bool);
		void SetEqType(const char*);
		
		/* This function will be used to call the model */
		Field(int,double*,double*,double*,double*,double*,double*,bool,bool);
		
	private:
		/* model parameters */
		double mui_,irho_,r0_,r1_,d_,xt_,xp_,dipshift_,diptilt_;
		const char *eqtype_;
		bool Edwards_;
		
		
		/* coordinate conversions for positions */
		void _SysIII2Mag(int,double*,double*,double*,double*,double*,double*,double*);
		void _PolSysIII2Mag(int,double*,double*,double*,double*,double*,double*,double*);
		
		
		/* coordinate conversion for magnetic field vector */
		
		/* Azimuthal field */
		_AzimuthalField(int,double*,double*);
		
};
