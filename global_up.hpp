/*
ver.03.03
file global_up.hpp
*/

#pragma once
#include <vector>
#include "common.h"
#include "global.h"

class Global_UP: virtual public Global{
public:
	Global_UP(
		const unsigned loadSteps,
		const double tolerance_equivalence,
		const unsigned max_iteration,
		const unsigned interval,
		const unsigned num_dof_p
	):
		Global(loadSteps,tolerance_equivalence,max_iteration,interval,false),
		num_dof_p(num_dof_p)
	{};
	virtual ~Global_UP(){};

	template<class PrePost_type,class Elem_type,class Material_type>
	void initialize(string problem,Material_type* material,const unsigned DIM,const unsigned NODE_IN_ELEM,const unsigned NUM_GP)
	{
		std::vector<Elem_type*> elem_buffer;
		prepost=new PrePost_type(
			vdof,num_dof_dirichlet,num_dof_nond,T,problem,elem_buffer,material,DIM,NODE_IN_ELEM,NUM_GP
		);
		initialize_unique(elem_buffer);
		after_preprocess();
		allocate_K();
	};

protected:
	template<class Elem_type>
	void initialize_unique(std::vector<Elem_type*>& elem_buffer){
		Element::Up_cast(Global::velements,elem_buffer);
		Element::Up_cast(velements,elem_buffer);

		int pressure_id=num_dof_nond;
		for(unsigned int i=0;i<velements.size();i++){
			velements[i]->put_pressureID(pressure_id);
			velements[i]->put_pressureID_(pressure_id+num_dof_dirichlet);
			pressure_id+=num_dof_p;
		}
	};

	void findK_unique(void){
		for(unsigned i=0;i<velements.size();i++)
			velements[i]->calcK_UP(*K);
	};

	void update_unique(void){
		for(unsigned i=0;i<velements.size();i++)
			velements[i]->calcT_p(T);
	};
	void affect_solution(double* delta_u){
		Global::affect_solution(delta_u);
		int pressure_id=num_dof_nond;
		for(int i=0;i<velements.size();i++)
			for(unsigned alpha=0;alpha<num_dof_p;alpha++)
				velements[i]->increment_pressure(alpha,delta_u[pressure_id++]);
	};

private:
	virtual void findK(void){
		using ABFEM::CONSTANTS::g_kpp_diag;
		Global::findK();
		findK_unique();
		for(unsigned i=num_dof_nond;i<num_row;i++)
			(*K)(i,i)+=g_kpp_diag;
#ifdef _DEBUG
		using namespace std;
		static bool isfirst=true;
		if(isfirst){
			cout << "printing K" << endl;
			isfirst=false;
			K->print("K.csv");
		}
#endif
	};

	virtual void update_iteration(void){
		Global::update_iteration();
		update_unique();
		F.ZeroOut();
		//issue12
		//for(unsigned int i=0;i<velements.size();i++)
			//velements[i]->calcBodyForce(F);
		
#ifdef _DEBUG
		static bool isfirst=true;
		if(isfirst){
			T.Print("T.csv");
			isfirst=false;
		}
#endif
	};

	void find_matrix_size(void){
		num_row_p=velements.size()*num_dof_p;
		num_row=num_dof_nond+num_row_p;
	};
	std::vector<Element_UP*> velements;
	unsigned num_row_p;
	const unsigned num_dof_p;
};

