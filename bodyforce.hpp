/*
ver.03.03
file bodyforce.hpp
*/

#pragma once
#include "global_pstab.hpp"

class Global_bodyforce: public Global_pstab{
public:
	Global_bodyforce(
		const unsigned loadSteps,
		const double tolerance_equivalence,
		const unsigned max_iteration,
		const unsigned interval,
		bool symmetric=false
	):
		Global_pstab(loadSteps,tolerance_equivalence,max_iteration,interval,symmetric),
		Global_UP(loadSteps,tolerance_equivalence,max_iteration,interval,1),
		Global(loadSteps,tolerance_equivalence,max_iteration,interval,symmetric)
	{
	};
	virtual ~Global_bodyforce(){};
protected:
	void update_iteration(void){
		F.ZeroOut();
		Global_pstab::update_iteration();
		for(unsigned int i=0;i<velements.size();i++)
			velements[i]->Tetra1::calcBodyForce(F);
	};
};

double g_density=0.000;
double g_gravity=9.81;
double g_bodyforce_normal[3]={0.0,0.0,-1.0};
int g_loadstep_cnt=1;
int g_loadsteps;
