#include "bessel.h"

void j0(int n, double *x, double *j) {
	
	/* define the parameters for this approx. */
	double lambda = 0.865;
	double q = 0.7172491568;
	double p0 = 0.6312725339;
	double pt0 = 0.4308049446;
	double p1 = 0.3500347951;
	double pt1 = 0.4678202347;
	double p2 = -0.06207747907;
	double pt2 = 0.04253832927;
	double l4 = lambda*lambda;
	l4 = l4*l4;
	
	/* some temporary variables */
	double x2, A, B, C, l4x2p1, rtl4x2p1, rtrtl4x2p1;
	
	/* loop through each point */
	int i;
	for (i=0;i<n;i++) {
		/* these temporary variables are used in multiple places so
		 * defining them once should hopefully speed things up a bit */
		x2 = x[i]*x[i];
		l4x2p1 = 1.0 + l4*x2;	
		rtl4x2p1 = sqrt(l4x2p1);
		rtrtl4x2p1 = sqrt(rtl4x2p1);
		
		A = 1.0/(rtrtl4x2p1*(1.0 + q*x2));
		B = p0 + p1*x2 + p2*rtl4x2p1;
		C = (pt0 + pt1*x2)*rtl4x2p1 + pt2*x2;
		
		j[i] = A*(B*cos(x[i]) + C*sin(x[i])/x[i]);
	}
		
}

void j0(int n, double *x, double multx, double *j) {
	
	/* define the parameters for this approx. */
	double lambda = 0.865;
	double q = 0.7172491568;
	double p0 = 0.6312725339;
	double pt0 = 0.4308049446;
	double p1 = 0.3500347951;
	double pt1 = 0.4678202347;
	double p2 = -0.06207747907;
	double pt2 = 0.04253832927;
	double l4 = lambda*lambda;
	l4 = l4*l4;
	
	/* some temporary variables */
	double mx, x2, A, B, C, l4x2p1, rtl4x2p1, rtrtl4x2p1;
	
	/* loop through each point */
	int i;
	for (i=0;i<n;i++) {
		/* these temporary variables are used in multiple places so
		 * defining them once should hopefully speed things up a bit */
		mx = x[i]*multx;
		x2 = mx*mx;
		l4x2p1 = 1.0 + l4*x2;	
		rtl4x2p1 = sqrt(l4x2p1);
		rtrtl4x2p1 = sqrt(rtl4x2p1);
		
		A = 1.0/(rtrtl4x2p1*(1.0 + q*x2));
		B = p0 + p1*x2 + p2*rtl4x2p1;
		C = (pt0 + pt1*x2)*rtl4x2p1 + pt2*x2;
		
		j[i] = A*(B*cos(mx) + C*sin(mx)/mx);
	}
		
}

void j1(int n, double *x, double *j) {
	
	/* parameters from the paper */
	double lambda = 0.1;
	double q1 = 0.4120981204;
	double q2 = 0.006571619275;
	double P0 = -0.7763224930;
	double p0 = 1.776322448;
	double P1 = -0.03147133771;
	double p1 = 0.2250803518;
	double P2 = -(2*pow(lambda,1.5))*q2/sqrt(M_PI);
	double p2 = (2*sqrt(lambda))*q2/sqrt(M_PI);
	double l2 = lambda*lambda;
	
	/* temporary variables */
	double x2, x4, A, B, C, ac, rtac, bc;
	
	/* loop through each value */
	int i;
	for (i=0;i<n;i++) {
		x2 = x[i]*x[i];
		x4 = x2*x2; 
		ac = sqrt(1 + l2*x2);
		rtac = sqrt(ac);
		A = 1.0/(2.0*rtac);
		bc = (1 + q1*x2 + q2*x4);
		B = (p0 + p1*x2 + p2*x4)/bc;
		C = (x[i]/ac)*(P0 + P1*x2 + P2*x4)/bc;
		j[i] = A*(B*sin(x[i]) + C*cos(x[i]));
	}

}

void j1(int n, double *x, double multx, double *j) {
	
	/* parameters from the paper */
	double lambda = 0.1;
	double q1 = 0.4120981204;
	double q2 = 0.006571619275;
	double P0 = -0.7763224930;
	double p0 = 1.776322448;
	double P1 = -0.03147133771;
	double p1 = 0.2250803518;
	double P2 = -(2*pow(lambda,1.5))*q2/sqrt(M_PI);
	double p2 = (2*sqrt(lambda))*q2/sqrt(M_PI);
	double l2 = lambda*lambda;
	
	/* temporary variables */
	double mx, x2, x4, A, B, C, ac, rtac, bc;
	
	/* loop through each value */
	int i;
	for (i=0;i<n;i++) {
		mx = multx*x[i];
		x2 = mx*mx;
		x4 = x2*x2; 
		ac = sqrt(1 + l2*x2);
		rtac = sqrt(ac);
		A = 1.0/(2.0*rtac);
		bc = (1 + q1*x2 + q2*x4);
		B = (p0 + p1*x2 + p2*x4)/bc;
		C = (mx/ac)*(P0 + P1*x2 + P2*x4)/bc;
		j[i] = A*(B*sin(mx) + C*cos(mx));
	}

}
