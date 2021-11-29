#include "internalmodel.h"

InternalModel::InternalModel() {
	
	/* load all of the models */
	vip4 = new Internal(&_binary_vip4coeffs_bin_start);
	jrm09 = new Internal(&_binary_jrm09coeffs_bin_start);
	gsfc13ev = new Internal(&_binary_gsfc13evcoeffs_bin_start);
	gsfc15ev = new Internal(&_binary_gsfc15evcoeffs_bin_start);
	gsfc15evs = new Internal(&_binary_gsfc15evscoeffs_bin_start);
	isaac = new Internal(&_binary_isaaccoeffs_bin_start);
	jpl15ev = new Internal(&_binary_jpl15evcoeffs_bin_start);
	jpl15evs = new Internal(&_binary_jpl15evscoeffs_bin_start);
	o4 = new Internal(&_binary_o4coeffs_bin_start);
	o6 = new Internal(&_binary_o6coeffs_bin_start);
	p11a = new Internal(&_binary_p11acoeffs_bin_start);
	sha = new Internal(&_binary_shacoeffs_bin_start);
	u17ev = new Internal(&_binary_u17evcoeffs_bin_start);
	v117ev = new Internal(&_binary_v117evcoeffs_bin_start);
	vipal = new Internal(&_binary_vipalcoeffs_bin_start);
	vit4 = new Internal(&_binary_vit4coeffs_bin_start);
		
	/* set the current model */
	CurrentModel_ = jrm09;
	strcpy(CurrentModelName_,"JRM09");
	
	/* default parameters */
	CartIn_ = true;
	CartOut_ = true;
	CurrentModel_->SetCartIn(CartIn_);
	CurrentModel_->SetCartOut(CartOut_);
	
}


InternalModel::~InternalModel() {
	
	/* delete all of the models */
	delete vip4;
	delete jrm09;
	delete gsfc13ev;
	delete gsfc15ev;
	delete gsfc15evs;
	delete isaac;
	delete jpl15ev;
	delete jpl15evs;
	delete o4;
	delete o6;
	delete p11a;
	delete sha;
	delete u17ev;
	delete v117ev;
	delete vipal;
	delete vit4;
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
		CurrentModel_ = jrm09;
	} else if (strcmp(ModelName,"VIP4") == 0) {
		CurrentModel_ = vip4;
	} else if (strcmp(ModelName,"GSFC13EV") == 0) {
		CurrentModel_ = gsfc13ev;
	} else if (strcmp(ModelName,"GSFC15EV") == 0) {
		CurrentModel_ = gsfc15ev;
	} else if (strcmp(ModelName,"GSFC15EVS") == 0) {
		CurrentModel_ = gsfc15evs;
	} else if (strcmp(ModelName,"ISAAC") == 0) {
		CurrentModel_ = isaac;
	} else if (strcmp(ModelName,"JPL15EV") == 0) {
		CurrentModel_ = jpl15ev;
	} else if (strcmp(ModelName,"JPL15EVS") == 0) {
		CurrentModel_ = jpl15evs;
	} else if (strcmp(ModelName,"O4") == 0) {
		CurrentModel_ = o4;
	} else if (strcmp(ModelName,"O6") == 0) {
		CurrentModel_ = o6;
	} else if (strcmp(ModelName,"P11A") == 0) {
		CurrentModel_ = p11a;
	} else if (strcmp(ModelName,"SHA") == 0) {
		CurrentModel_ = sha;
	} else if (strcmp(ModelName,"U17EV") == 0) {
		CurrentModel_ = u17ev;
	} else if (strcmp(ModelName,"V117EV") == 0) {
		CurrentModel_ = v117ev;
	} else if (strcmp(ModelName,"VIPAL") == 0) {
		CurrentModel_ = vipal;
	} else if (strcmp(ModelName,"VIT4") == 0) {
		CurrentModel_ = vit4;
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

void InternalModel::Field(	double p0, double p1, double p2,
							double *B0, double *B1, double *B2) {
	
	CurrentModel_->Field(p0,p1,p2,B0,B1,B2);
								
}
