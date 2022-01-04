#include "internalmodel.h"

InternalModel::InternalModel() {
	
	/* load all of the models */
	vip4_ = &vip4;
	jrm09_ = &jrm09;
	jrm33_ = &jrm33;
	gsfc13ev_ = &gsfc13ev;
	gsfc15ev_ = &gsfc15ev;
	gsfc15evs_ = &gsfc15evs;
	isaac_ = &isaac;
	jpl15ev_ = &jpl15ev;
	jpl15evs_ = &jpl15evs;
	o4_ = &o4;
	o6_ = &o6;
	p11a_ = &p11a;
	sha_ = &sha;
	u17ev_ = &u17ev;
	v117ev_ = &v117ev;
	vipal_ = &vipal;
	vit4_ = &vit4;
		
	/* set the current model */
	CurrentModel_ = jrm09_;
	strcpy(CurrentModelName_,"JRM09");
	
	/* default parameters */
	CartIn_ = true;
	CartOut_ = true;
	CurrentModel_->SetCartIn(CartIn_);
	CurrentModel_->SetCartOut(CartOut_);
	
}


InternalModel::~InternalModel() {

}	

void InternalModel::SetCartIn(bool CartIn){
	
	CartIn_ = CartIn;
	CurrentModel_->SetCartIn(CartIn_);
}

bool InternalModel::GetCartIn() {
	
	return CartIn_;
}

void InternalModel::SetCartOut(bool CartOut){
	
	CartOut_ = CartOut;
	CurrentModel_->SetCartOut(CartOut_);
}

bool InternalModel::GetCartOut() {
	
	return CartOut_;
}

void InternalModel::SetModel(char *ModelName) {

	/* Find the correct model and set it */
	bool validmodel = true;
	if (strcmp(ModelName,"JRM09") == 0) {
		CurrentModel_ = jrm09_;
	} else if (strcmp(ModelName,"JRM33") == 0) {
		CurrentModel_ = jrm33_;
	} else if (strcmp(ModelName,"VIP4") == 0) {
		CurrentModel_ = vip4_;
	} else if (strcmp(ModelName,"GSFC13EV") == 0) {
		CurrentModel_ = gsfc13ev_;
	} else if (strcmp(ModelName,"GSFC15EV") == 0) {
		CurrentModel_ = gsfc15ev_;
	} else if (strcmp(ModelName,"GSFC15EVS") == 0) {
		CurrentModel_ = gsfc15evs_;
	} else if (strcmp(ModelName,"ISAAC") == 0) {
		CurrentModel_ = isaac_;
	} else if (strcmp(ModelName,"JPL15EV") == 0) {
		CurrentModel_ = jpl15ev_;
	} else if (strcmp(ModelName,"JPL15EVS") == 0) {
		CurrentModel_ = jpl15evs_;
	} else if (strcmp(ModelName,"O4") == 0) {
		CurrentModel_ = o4_;
	} else if (strcmp(ModelName,"O6") == 0) {
		CurrentModel_ = o6_;
	} else if (strcmp(ModelName,"P11A") == 0) {
		CurrentModel_ = p11a_;
	} else if (strcmp(ModelName,"SHA") == 0) {
		CurrentModel_ = sha_;
	} else if (strcmp(ModelName,"U17EV") == 0) {
		CurrentModel_ = u17ev_;
	} else if (strcmp(ModelName,"V117EV") == 0) {
		CurrentModel_ = v117ev_;
	} else if (strcmp(ModelName,"VIPAL") == 0) {
		CurrentModel_ = vipal_;
	} else if (strcmp(ModelName,"VIT4") == 0) {
		CurrentModel_ = vit4_;
	} else {
		printf("Invalid model name: %s, ignoring...\n",ModelName);
		validmodel = false;
	}
	
	/* change the stored model name and set model parameters */
	if (validmodel) {
		strcpy(CurrentModelName_,ModelName);
		CurrentModel_->SetCartIn(CartIn_);
		CurrentModel_->SetCartOut(CartOut_);
	}
}

void InternalModel::GetModel(char *ModelName) {
	
	strcpy(ModelName,CurrentModelName_);
}

void InternalModel::Field(int n, double *p0, double *p1, double *p2,
							double *B0, double *B1, double *B2) {
	
	CurrentModel_->Field(n,p0,p1,p2,B0,B1,B2);
								
}

void InternalModel::Field(int n, double *p0, double *p1, double *p2,
							int MaxDeg, double *B0, double *B1, double *B2) {
	
	CurrentModel_->Field(n,p0,p1,p2,MaxDeg,B0,B1,B2);
								
}

void InternalModel::Field(	double p0, double p1, double p2,
							double *B0, double *B1, double *B2) {
	
	CurrentModel_->Field(p0,p1,p2,B0,B1,B2);
								
}

void InternalModel::Field(	double p0, double p1, double p2, int MaxDeg,
							double *B0, double *B1, double *B2) {
	
	CurrentModel_->Field(p0,p1,p2,MaxDeg,B0,B1,B2);
								
}
