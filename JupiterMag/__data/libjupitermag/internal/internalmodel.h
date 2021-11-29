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
		void Field(int,double*,double*,double*,double*,double*,double*);
		void Field(double,double,double,double*,double*,double*);
				
		/* these objects are the models to use */
		Internal *vip4, *jrm09, *gsfc13ev, *gsfc15ev, *gsfc15evs, *isaac,
				*jpl15ev, *jpl15evs, *o4, *o6, *p11a, *sha, *u17ev, *v117ev,
				*vipal, *vit4;
	private:
		Internal *CurrentModel_;
		char CurrentModelName_[32];
		bool CartIn_;
		bool CartOut_;
};


#endif
