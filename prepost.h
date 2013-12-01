/*
ver.03.03
file prepost.h
*/
#pragma once
#include <sstream>
#include <string>
#include <vector>
#include "element.h"

class PrePostProcessor{
public:
	PrePostProcessor(
		std::vector<Dof*>& vdof,
		unsigned& num_dof_dirichlet,
		unsigned& num_dof_nond,
		MATRIX& T,
		std::string projectname,
		const unsigned DIM,
		const unsigned NODE_IN_ELEM,
		const unsigned NUM_GP
	):
		//num_nodes(global.num_nodes),
		vdof(vdof),
		num_dof_dirichlet(num_dof_dirichlet),
		num_dof_nond(num_dof_nond),
		T(T),
		projectname(projectname),
		DIM(DIM),
		NODE_IN_ELEM(NODE_IN_ELEM),
		NUM_GP(NUM_GP)
	{
		m_output_element << "gid-MAFINS/"
			<< projectname << "-element.gid/"
			<< projectname << "-element.post.res";
	};
	virtual void post_process(unsigned loadstep)=0;
	virtual void save_result(unsigned loadstep)=0;
	virtual ~PrePostProcessor(){};

	void debug(void){
		using namespace std;
		for(unsigned int i=0;i<vdof.size()/3;i++){
			cout << i << "\t";
			cout <<
				vdof[raw_to_mi[i*3  ]]->X0 << "\t" << 
				vdof[raw_to_mi[i*3+1]]->X0 << "\t" << 
				vdof[raw_to_mi[i*3+2]]->X0 << "\t";
			cout <<
				vdof[raw_to_mi[i*3  ]]->X << "\t" << 
				vdof[raw_to_mi[i*3+1]]->X << "\t" << 
				vdof[raw_to_mi[i*3+2]]->X <<endl;
		}
	};

protected:
	std::vector<Element_output*> velements;
	//Node*& raw_nodes;
	Node* raw_nodes;
	int* raw_to_mi;
	unsigned num_nodes;

	//from global
	std::vector<Dof*>& vdof;
	unsigned& num_dof_dirichlet;
	unsigned& num_dof_nond;
	MATRIX& T;

	const unsigned DIM;
	const unsigned NODE_IN_ELEM;
	const unsigned NUM_GP;
	std::string projectname;
	std::ostringstream m_output_element;
};

