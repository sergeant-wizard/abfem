/*
ver.03.03
file global_pstab.hpp
*/

#pragma once
#include "global_up.hpp"
#include "elemface.h"

class Global_pstab: virtual public Global_UP{
public:
	Global_pstab(
		const unsigned loadSteps,
		const double tolerance_equivalence,
		const unsigned max_iteration,
		const unsigned interval,
		bool symmetric=false
	):
		Global_UP(loadSteps,tolerance_equivalence,max_iteration,interval,1),
		Global(loadSteps,tolerance_equivalence,max_iteration,interval,symmetric)
	{
		for(unsigned i=0;i<velem_faces.size();i++)
			delete velem_faces[i];
	};
	virtual ~Global_pstab(){};

	template<class PrePost_type,class Material_type>
	void initialize(string problem,Material_type* material,string elemface_file){
		prepost=new PrePost_type(
			vdof,num_dof_dirichlet,num_dof_nond,T,problem,velements,material,3,4,1
		);
		Global_UP::initialize_unique<Tetra1_UP>(velements);
		after_preprocess();

		std::vector<Elem_face*> elemface_buffer;
		initialize_unique<Elem_face>(velements,elemface_file,elemface_buffer);

		allocate_K();
	};

	void update_iteration(void){
		Global::update_iteration();
		Global_UP::update_unique();
		for(unsigned int i=0;i<velem_faces.size();i++)
			velem_faces[i]->update_iteration();
		for(unsigned int i=0;i<velem_faces.size();i++)
			velem_faces[i]->calcT(T);
	};

	virtual void findK(void){
		Global::findK();
		Global_UP::findK_unique();
		for(unsigned int i=0;i<velem_faces.size();i++)
			velem_faces[i]->calcK(*K);

	};
protected:
	template<class ElemFace_type,class Element_type>
	void initialize_unique(
		std::vector<Element_type*>& elem_buffer,
		string elemface_file,
		std::vector<ElemFace_type*>& velem_faces_buffer
	){
		static const int faces[4][3]={{1,2,3},{0,2,3},{0,1,3},{0,1,2}};
		for(unsigned i=0;i<elem_buffer.size();i++){
			for(int j=0;j<4;j++){
				ElemFace_type* newef=new ElemFace_type(
					elem_buffer[i]->get_node(faces[j][0]),
					elem_buffer[i]->get_node(faces[j][1]),
					elem_buffer[i]->get_node(faces[j][2]),
					elem_buffer[i],
					elemface_file
				);
				bool flag=false;
				for(unsigned k=0;k<velem_faces_buffer.size();k++)
					if(*velem_faces_buffer[k] == *newef){
						*velem_faces_buffer[k]+=*newef;
						velem_faces_buffer.push_back(velem_faces_buffer[k]);
						velem_faces_buffer.erase(velem_faces_buffer.begin()+k);
						flag=true;
						break;
					}
				if(!flag)
					velem_faces_buffer.insert(velem_faces_buffer.begin(),newef);
				else
					delete newef;
			}
		}
		Element::Up_cast(velem_faces,velem_faces_buffer);
		for(unsigned i=0;i<velem_faces.size();i++)
			velem_faces[i]->initialize(*K);
	};

	std::vector<Tetra1_UP*> velements;
	std::vector<Elem_face*> velem_faces;
private:
};

