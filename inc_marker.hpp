/*
ver.03.03
file inc_marker.hpp
*/

#pragma once
#include "tetra1_marker.h"
#include "tetra1_up.hpp"
class Inc_marker_element:public Tetra1_UP,public Tetra1_marker{
public:
	Inc_marker_element(
		int* node_index,
		Node* raw_nodes,
		Default_material* material
	):
		Tetra1(node_index,raw_nodes,material),
		Tetra1_UP(node_index,raw_nodes,material),
		Tetra1_marker(node_index,raw_nodes,material),
		Element_UP(1)
	{};
	void initialize_sparse(Sparse& K){
		Tetra1_UP::initialize_sparse(K);
	};
	void update_loadstep(void){
		Tetra1_marker::update_loadstep();
		m_p[0]=0;
	};
	void update_iteration(void){
		Tetra1::update_F();
		Tetra1::update_dNdx();
		material->getC_e(m_c_hat,&m_F_rel_);
		material->getSigma(&m_sigma_hat,&m_F_rel_);
		Tetra1_UP::update_iteration_unique();

		update_markers();
		m_marker_volume=m_marker_Volume*m_volume/m_Volume;
	};
	void calcK(Sparse& K){
		Tetra1_marker::calcK(K);
	};
	void calcT(MATRIX& T){
		Tetra1_marker::calcT(T);
	};

};

#include "elemface.h"

class Elem_face_marker:public Elem_face{
	typedef Elem_face Parent;
	typedef Tetra1_marker INPUT_ELEMENT;
public:
	Elem_face_marker(Node* a,Node* b,Node* c,Inc_marker_element* elem,string filename):
		Parent(a,b,c,elem,filename),
		elem0(elem),
		elem1(NULL)
	{};
	void operator+=(Elem_face_marker& right){
		Parent::elem1=right.Parent::elem0;
		elem1=right.elem0;
	};
	void update_iteration(void){
		//using ABFEM::CONSTANTS::g_delta;
		//using ABFEM::CONSTANTS::g_delta_weak;
		if(!isInterior())
			return;
		Parent::update_iteration();

#ifdef BETA
		//issue11
		if(elem0->effective_flag && elem1->effective_flag){
			if(elem0->isCut() && elem1->isCut())
				coef=area*m_delta_h_weak;
			else if(
				(elem0->isCut() && elem1->get_weight()>0) ||
				(elem1->isCut() && elem0->get_weight()>0)
			)
				coef=0;
			else
				coef=area*m_delta_h;
		}else
			coef=0;
		return;
#endif

		if(elem0->isCut() && elem1->isCut()){
			coef=area*m_delta_h_weak;
			//coef=area*m_delta_h;
		}else if(
			//border & outside
			(
				(elem0->get_weight()==0 && elem1->get_weight()>0)||
				(elem1->get_weight()==0 && elem0->get_weight()>0)
			)||(
			//border & inside
				(elem0->isCut() && elem1->get_weight()>0) ||
				(elem1->isCut() && elem0->get_weight()>0)
			)
		){
			coef=0;
		}else{
			coef=area*m_delta_h;
		}
	};

	Tetra1_marker* get_elem(unsigned elem_index){
		if(!elem_index)
			return elem0;
		else
			return elem1;
	};

private:
	Tetra1_marker* elem0;
	Tetra1_marker* elem1;
};

#include "global_pstab.hpp"
#include "global_marker.hpp"

class Global_inc_marker:virtual public Global_pstab,virtual public Global_marker{
public:
	Global_inc_marker(
		const unsigned loadSteps,
		const double tolerance_equivalence,
		const unsigned max_iteration,
		const unsigned interval,
		const double dirichlet,
		bool symmetric=false
	):
		Global(loadSteps,tolerance_equivalence,max_iteration,interval),
		Global_UP(loadSteps,tolerance_equivalence,max_iteration,interval,1),
		Global_pstab(loadSteps,tolerance_equivalence,max_iteration,interval,dirichlet),
		Global_marker(loadSteps,tolerance_equivalence,max_iteration,interval,dirichlet)
	{};

	template<class PrePost_type,class Material_type>
	void initialize(
		string projectname,Material_type* front_material,
		const unsigned DIM,const unsigned NODE_IN_ELEM,
		Material_type* back_material,
		string elemface_file
	){
		std::vector<Inc_marker_element*> elem_buffer;
		prepost=new PrePost_type(
			vdof,num_dof_dirichlet,num_dof_nond,T,projectname,elem_buffer,
			front_material,3,4,1,markers_handler,back_material
		);

		Global_UP::initialize_unique(elem_buffer);
		Element::Up_cast(Global_marker::velements,elem_buffer);

		after_preprocess();
		std::vector<Elem_face_marker*> velem_faces_buffer;
		Global_pstab::initialize_unique(elem_buffer,elemface_file,velem_faces_buffer);
		allocate_K();
	};
};

