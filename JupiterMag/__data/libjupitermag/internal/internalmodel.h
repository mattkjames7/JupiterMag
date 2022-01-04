#ifndef __INTERNALMODEL_H__
#define __INTERNALMODEL_H__
#include <stdio.h>
#include <stdlib.h>
#include "internal.h"
#include <string.h>

class InternalModel {
	
	public:
		/* constructor */
		InternalModel();
		
		/* destructor */
		~InternalModel();
		
		/* set model parameters */
		void SetCartIn(bool);
		void SetCartOut(bool);
		bool GetCartIn();
		bool GetCartOut();
		void SetModel(char *);
		void GetModel(char *);

		/* Field functions */
		void Field(int,double*,double*,double*,int,double*,double*,double*);
		void Field(int,double*,double*,double*,double*,double*,double*);
		void Field(double,double,double,int,double*,double*,double*);
		void Field(double,double,double,double*,double*,double*);
				
		/* these objects are the models to use */
		Internal *vip4_, *jrm09_, *jrm33_, *gsfc13ev_, *gsfc15ev_, *gsfc15evs_, 
				*isaac_, *jpl15ev_, *jpl15evs_, *o4_, *o6_, *p11a_, 
				*sha_, *u17ev_, *v117ev_, *vipal_, *vit4_;
	private:
		Internal *CurrentModel_;
		char CurrentModelName_[32];
		bool CartIn_;
		bool CartOut_;
};


#endif
