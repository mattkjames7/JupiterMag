#include "libinternal.h"

void InternalField(int l, double *p0, double *p1, double *p2,
					double *B0, double *B1, double *B2,
					bool PolIn, bool PolOut, const char *Model) {
	
	/* select the appropriate model (default to VIP4) */
	if (strcmp(Model,"VIP4") == 0) {
		vip4.Field(l,p0,p1,p2,B0,B1,B2,PolIn,PolOut);
	} else if (strcmp(Model,"JRM09") == 0) {
		jrm09.Field(l,p0,p1,p2,B0,B1,B2,PolIn,PolOut);
	} else {
		vip4.Field(l,p0,p1,p2,B0,B1,B2,PolIn,PolOut);
	}		
						
}
