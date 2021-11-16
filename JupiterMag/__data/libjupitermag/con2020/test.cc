#include "test.h"

int main() {
	
	/* set model to integral */
	con2020.SetEqType("integral");
	
	/* call the model */
	double x = 5.0;
	double y = 0.0;
	double z = 0.0;
	double Bx, By, Bz;
	con2020.Field(1,&x,&y,&z,&Bx,&By,&Bz);
	
	printf("B = [%f, %f, %f]\n",Bx,By,Bz);
	
}
