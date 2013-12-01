/*
ver.03.03
file gid.cpp
*/

#include "gid.h"
#include "gid.hpp"

using std::endl;

const std::string Gid::PATH="gid/";

Gid::~Gid(){
	for(unsigned int i=0;i<velements.size();i++)
		delete velements[i];
	velements.clear();
	safe_delete(raw_nodes);
	safe_delete(raw_to_mi);
}

void Gid::save_result(unsigned loadstep){
	m_result.str("");

	gid_result_<Dof>      (m_result,loadstep,"displacement"  ,&Dof::get_displacement);
	gid_result<MATRIX>   (m_result,loadstep,"internal force",T);
	gid_result           (m_result,loadstep,"stress_xx"     ,&Element_output::getStress,0,0);
	gid_result           (m_result,loadstep,"stress_yy"     ,&Element_output::getStress,1,1);
	gid_result           (m_result,loadstep,"stress_xy"     ,&Element_output::getStress,0,1);
	if(DIM==3){
		gid_result(m_result,loadstep,"stress_zz",&Element_output::getStress,2,2);
		gid_result(m_result,loadstep,"stress_xz",&Element_output::getStress,0,2);
		gid_result(m_result,loadstep,"stress_yz",&Element_output::getStress,1,2);
	}
	gid_result(m_result,loadstep,"deform_xx"     ,&Element_output::get_F_,0,0);
	gid_result(m_result,loadstep,"deform_xy"     ,&Element_output::get_F_,0,1);
	gid_result(m_result,loadstep,"deform_yx"     ,&Element_output::get_F_,1,0);
	gid_result(m_result,loadstep,"deform_yy"     ,&Element_output::get_F_,1,1);
	if(DIM==3){
		gid_result(m_result,loadstep,"deform_xz",&Element_output::get_F_,0,2);
		gid_result(m_result,loadstep,"deform_yz",&Element_output::get_F_,1,2);
		gid_result(m_result,loadstep,"deform_zx",&Element_output::get_F_,2,0);
		gid_result(m_result,loadstep,"deform_zy",&Element_output::get_F_,2,1);
		gid_result(m_result,loadstep,"deform_zz",&Element_output::get_F_,2,2);
	}
	gid_result(m_result,loadstep,"volumetric strain",&Element_output::get_J_);
	gid_result(m_result,loadstep,"pressure",&Element_output::get_pressure);
	gid_result(m_result,loadstep,"potential",&Element_output::get_potential);
	gid_result(m_result,loadstep,"volume",&Element_output::get_volume);

	//issue11
	//gid_result(m_result,loadstep,"effective_flag",&Element_output::get_effective_flag);
}

void Gid::post_process(unsigned loadstep){
	static bool isfirst=true;
	ofstream fout;
	if(isfirst){
		isfirst=false;
		fout.open(m_output_element.str().c_str());
		fout << "GiD Post Results File 1.0" <<endl<<endl;
		fout << "GaussPoints \"gauss_points\" ElemType " << gid_elemtype << endl;
		fout << "Number of Gauss Points: " << NUM_GP << endl;
		fout << "Nodes not included" << endl;
		fout << "Natural Coordinates: Internal" << endl;
		fout << "End Gausspoints" << endl;
	}else{
		fout.open(m_output_element.str().c_str(),std::ios::app);
	}

	fout << m_result.str();

	/*
	//residue
	fout << "Result \"residue force\" \"FEM\" ";
	fout << loadstep << " Vector OnNodes" << endl;
	fout << "Values" <<endl;
	for(int i=0;i<num_nodes;i++){
		fout << (i+1) << "\t";
		for(int j=0;j<3;j++){
			MATRIXINDEX mi=raw_to_mi[i*DIM+j];
			if(mi<num_dof_dirichlet)
				fout << "0 ";
			else
				fout << -Rminus(mi-num_dof_dirichlet,0) << " ";
		}
		fout << endl;
	}
	fout << "End Values" <<endl;
	*/
}

void Gid::gid_result(
	ostream& fout,unsigned loadstep,string name,
	double (Element_output::*func)(unsigned,unsigned,unsigned)const,unsigned i_,unsigned j_
){
	fout << "Result \"" << name << "\" \"FEM\" ";
	fout << loadstep << " ";
	fout << "Scalar OnGaussPoints gauss_points" << endl;
	fout << "Values" <<endl;
	for(unsigned i=0;i<velements.size();i++){
		fout << (i+1) << "\t";
		for(unsigned gp=0;gp<NUM_GP;gp++)
			fout << (velements[i]->*func)(gp,i_,j_) << endl;
	}
	fout << "End Values" <<endl;
}

void Gid::get_type(void){
	using namespace std;
	switch(DIM){
		default:
			switch(NODE_IN_ELEM){
				case 27:
					gid_elemtype="Hexahedra";
					break;
				case 8:
					gid_elemtype="Tetrahedra";
					break;
				default:
					gid_elemtype="Tetrahedra";
					break;
			}
			break;
	}
}

