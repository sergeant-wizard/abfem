/*
ver.03.03
file global.h
*/

#pragma once

#include <iomanip>
#include <iostream>
#include <string>
#include "common.h"
#include "matrix_wiz.h"
#include "sparse.h"
#include "element.h"
#include "prepost.h"

class Global{
public:
	typedef std::string string;
	typedef MATRIX_WIZ Default_matrix;
	typedef Sparse_wiz Default_sparse;
	Global(
		const unsigned loadSteps,
		const double tolerance_equivalence,
		const unsigned max_iteration,
		const unsigned interval,
		bool symmetric=true
	);
	virtual ~Global();
	bool solve(void);
	
	template<class PrePost_type,class Elem_type,class Material_type>
	void initialize(string projectname,Material_type* material,const unsigned DIM,const unsigned NODE_IN_ELEM,const unsigned NUM_GP);

	double& get_T(int row){
		return T(row+num_dof_dirichlet,0);
	};
protected:

	virtual void after_preprocess(void);
	virtual void update_loadstep(void);
	virtual void update_iteration(void);
	virtual void after_iteration(unsigned load_step);
	virtual void allocate_K(void);
	virtual void findK(void);

	virtual void find_matrix_size(void){
		num_row=num_dof_nond;
	};
	virtual void affect_solution(double* delta_u){
		for(int i=0;i<num_dof_nond;i++){
			vdof[i+num_dof_dirichlet]->x+=delta_u[i];
			vdof[i+num_dof_dirichlet]->u+=delta_u[i];
		}
	};

	//data(reference)
	std::vector<Dof*> vdof;
	std::vector<Element*> velements;

	PrePostProcessor* prepost;
	Default_sparse* K;
	Default_matrix T;
	Default_matrix F;
	Default_matrix Rminus;

	//unsigned num_nodes;
	unsigned num_dof_dirichlet;
	unsigned num_dof_nond;
	unsigned num_row;

	const unsigned loadSteps;
	const unsigned interval;
private:
	void iteration(void);
	void findT(void);
	void calculate(void);

	//options
	const unsigned max_iteration;
	const double tolerance_equivalence;
	const bool symmetric;
};

template<class PrePost_type,class Elem_type,class Material_type>
void Global::initialize(string problem,Material_type* material,const unsigned DIM,const unsigned NODE_IN_ELEM,const unsigned NUM_GP){
	std::vector<Elem_type*> elem_buffer;
	prepost=new PrePost_type(
		vdof,num_dof_dirichlet,num_dof_nond,T,problem,elem_buffer,material,DIM,NODE_IN_ELEM,NUM_GP
	);
	Element::Up_cast(velements,elem_buffer);

	after_preprocess();
	allocate_K();
}

