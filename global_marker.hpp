/*
ver.03.03
file global_marker.hpp
*/

#pragma once
#include "global.h"
#include "global_up.hpp"
#include "marker.h"
#include "tetra1_marker.h"

class Global_marker: virtual public Global{
public:
	Global_marker(
		const unsigned loadSteps,
		const double tolerance_equivalence,
		const unsigned max_iteration,
		const unsigned interval,
		const double dirichlet,
		bool symmetric=true
	):
		Global(loadSteps,tolerance_equivalence,max_iteration,interval,symmetric),
		markers_handler(NULL),
		dirichlet(dirichlet)
	{
	};
	virtual ~Global_marker(){
	};

	template<class PrePost_type,class Elem_type,class Material_type>
	void initialize(
		string projectname,Material_type* front_material,
		const unsigned DIM,const unsigned NODE_IN_ELEM,const unsigned NUM_GP,
		Material_type* back_material
	){
		std::vector<Elem_type*> elem_buffer;
		prepost=new PrePost_type(
			vdof,num_dof_dirichlet,num_dof_nond,T,projectname,elem_buffer,
			front_material,DIM,NODE_IN_ELEM,NUM_GP,markers_handler,back_material
		);
		Element::Up_cast(Global::velements,elem_buffer);
		Element::Up_cast(velements,elem_buffer);

		after_preprocess();
		allocate_K();
	};

protected:
	Markers_handler* markers_handler;
	std::vector<Element_marker*> velements;
	void update_loadstep(void){
		using std::cout;
		for(unsigned i=0;i<velements.size();i++)
			velements[i]->update_loadstep();
		cout << "handling markers...";
		cout.flush();
		markers_handler->update_element();
		for(unsigned int i=0;i<velements.size();i++)
			velements[i]->calcVolume();
		cout << "\r\x1b[K";
		Global::update_loadstep();
	};
private:

	void after_iteration(unsigned load_step){
		for(unsigned int i=0;i<velements.size();i++)
			velements[i]->move_markers();
		Global::after_iteration(load_step);
		meshcontrol();
	};

	void meshcontrol(void){
		/*
		//for shear
		static const double a_z=0.0;
		static const double b_z=1.0;
		static int cnt=1;
		for(int i=0;i<num_dof_dirichlet+num_dof_nond;i++){
			vdof[i]->x=vdof[i]->X=vdof[i]->X0;
			vdof[i]->u=0;
		}
		cnt++;
		*/
		//for cuffs
		static const double a_z=0.0;
		static const double b_z=1.0;
		static const double delta=dirichlet/(double)loadSteps;
		static int cnt=1;
		for(int i=0;i<num_dof_dirichlet+num_dof_nond;i++){
			if(vdof[i]->rawindex%3==2){
				if(vdof[i]->X0 <= b_z){
					vdof[i]->X=
						vdof[i]->X0+
						(vdof[i]->X0-a_z)/(b_z-a_z)*delta*(double)cnt;
				}else{
					vdof[i]->X=
						vdof[i]->X0+
						delta*(double)cnt;
				}
			}
			vdof[i]->x=vdof[i]->X;
			vdof[i]->u=0;
		}
		cnt++;
		/*
		//for ppm
		static const double a_z=0.0;
		static const double b_z=1.0;
		static const double c_z=1.5;
		static const double delta=dirichlet/(double)loadSteps;
		static int cnt=1;
		for(int i=0;i<num_dof_dirichlet+num_dof_nond;i++){
			if(vdof[i]->rawindex%3==2){
				if(vdof[i]->X0 <= b_z){
					vdof[i]->X=
						vdof[i]->X0+
						(vdof[i]->X0-a_z)/(b_z-a_z)*delta*(double)cnt;
				}else{
					vdof[i]->X=
						vdof[i]->X0+
						(c_z-vdof[i]->X0)/(c_z-b_z)*delta*(double)cnt;
				}
			}
			vdof[i]->x=vdof[i]->X;
			vdof[i]->u=0;
		}
		cnt++;
		*/
		/*ua series
		static const double min_z=0;
		static const double max_z=1;
		static const double delta=dirichlet/(double)loadSteps;
		static int cnt=1;
		for(int i=0;i<num_dof_dirichlet+num_dof_nond;i++){
			if(vdof[i]->rawindex%3==2){
				vdof[i]->X=
					vdof[i]->X0+
					(vdof[i]->X0-min_z)/(max_z-min_z)*delta*(double)cnt;
			}
			vdof[i]->x=vdof[i]->X;
			vdof[i]->u=0;
		}
		cnt++;
		*/
	};

	const double dirichlet;
};

