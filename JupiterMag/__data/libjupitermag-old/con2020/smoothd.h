#ifndef __SMOOTHD_H__
#define __SMOOTHD_H__
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/***********************************************************************
 * NAME : smoothd(z,dz,d)
 * 
 * DESCRIPTION : Smooth fucntion for crossing the current sheet 
 * (replaces the last bit of equation 12 in Edwards et al 2000).
 * 
 * INPUTS : 
 * 		double z	z-coordinate in dipole coordinate system (Rj)
 * 		double dz	Scale of the transition to use (Rj)
 * 		double d	Half thickness of the current sheet.
 * 
 * RETURNS : 
 * 		double out	Smoothed function across the current sheet.
 * 
 * ********************************************************************/
double smoothd(double z, double dz, double d);

#endif
