/*
ver.03.03
file gid_marker.hpp
*/

#include "gid.h"
#include "marker.h"

class Gid_marker:public Gid{
public:
	template<class Elem_type,class Material_type>
	Gid_marker(
		std::vector<Dof*>& vdof,
		unsigned& num_dof_dirichlet,
		unsigned& num_dof_nond,
		MATRIX& T,
		string projectname,
		std::vector<Elem_type*>& elem_buffer,
		Material_type* front_material,
		const unsigned DIM,
		const unsigned NODE_IN_ELEM,
		const unsigned NUM_GP,
		Markers_handler*& markers_handler,
		Material_type* back_material
	):
		PrePostProcessor(vdof,num_dof_dirichlet,num_dof_nond,T,projectname,DIM,NODE_IN_ELEM,NUM_GP),
		Gid(vdof,num_dof_dirichlet,num_dof_nond,T,projectname,elem_buffer,back_material,DIM,NODE_IN_ELEM,NUM_GP),
		markers_handler(markers_handler)
	{
		ostringstream marker_file;
		marker_file <<
			PATH <<
			projectname << "-marker.gid/" <<
			projectname << "-marker.dat";
		marker_file_out <<
			PATH << 
			projectname << "-marker.gid/" <<
			projectname << "-marker.post.res";
		Element::Up_cast(velements,elem_buffer);
		markers_handler=new Markers_handler(
			marker_file,
			elem_buffer,
			front_material
		);
	};
	virtual ~Gid_marker(){
		if(markers_handler)
			delete markers_handler;
	};

	void save_result(unsigned loadstep){
		m_result.str("");
		//parent
		Gid::save_result(loadstep);
		Gid::gid_result<Element_marker>(Gid::m_result,loadstep,"marker weight",&Element_marker::get_weight,velements);
		Gid::gid_result<Element_marker>(Gid::m_result,loadstep,"num marker",&Element_marker::get_num_marker,velements);

		gid_result(m_result,loadstep,"displacement",&Marker::getDisplacement);
		gid_result(m_result,loadstep,"stress_xx",&Marker::getStress,0,0);
		gid_result(m_result,loadstep,"stress_xy",&Marker::getStress,0,1);
		gid_result(m_result,loadstep,"stress_xz",&Marker::getStress,0,2);
		gid_result(m_result,loadstep,"stress_yy",&Marker::getStress,1,1);
		gid_result(m_result,loadstep,"stress_yz",&Marker::getStress,1,2);
		gid_result(m_result,loadstep,"stress_zz",&Marker::getStress,2,2);
		gid_result(m_result,loadstep,"deform_xx",&Marker::getStrain,0,0);
		gid_result(m_result,loadstep,"deform_xy",&Marker::getStrain,0,1);
		gid_result(m_result,loadstep,"deform_xz",&Marker::getStrain,0,2);
		gid_result(m_result,loadstep,"deform_yx",&Marker::getStrain,1,0);
		gid_result(m_result,loadstep,"deform_yy",&Marker::getStrain,1,1);
		gid_result(m_result,loadstep,"deform_yz",&Marker::getStrain,1,2);
		gid_result(m_result,loadstep,"deform_zx",&Marker::getStrain,2,0);
		gid_result(m_result,loadstep,"deform_zy",&Marker::getStrain,2,1);
		gid_result(m_result,loadstep,"deform_zz",&Marker::getStrain,2,2);
		gid_result<int>(m_result,loadstep,"parent_elem",&Marker::get_parent);
		gid_result<double>(m_result,loadstep,"potential",&Marker::getPotential);
	};

	void post_process(unsigned loadstep){
		using std::endl;
		Gid::post_process(loadstep);
		
		ofstream fout;
		static bool isfirst=true;
		if(isfirst){
			isfirst=false;
			fout.open(marker_file_out.str().c_str());
			fout << "GiD Post Results File 1.0" <<endl<<endl;
			fout << "GaussPoints \"gauss_points\" ElemType " << gid_elemtype << endl;
			fout << "Number of Gauss Points: 1" << endl;
			fout << "Nodes not included" << endl;
			fout << "Natural Coordinates: Internal" << endl;
			fout << "End Gausspoints" << endl;
		}else{
			fout.open(marker_file_out.str().c_str(),std::ios::app);
		}
		fout << m_result.str();
	};
private:
	void gid_result(ostream& fout,unsigned loadstep,string name,double (Marker::*func)(unsigned)const){
		using std::endl;
		fout << "Result \"" << name << "\" \"FEM\" ";
		fout << loadstep << " ";
		fout << "Vector OnNodes" << endl;
		fout << "Values" <<endl;
		for(unsigned i=0;i<markers_handler->num_markers;i++){
			fout << (i+1) << "\t";
			for(unsigned j=0;j<DIM;j++){
				fout << ((*markers_handler)(i)->*func)(j) << " ";
			}
			fout<<endl;
		}
		fout << "End Values" <<endl;
	};
	void gid_result(ostream& fout,unsigned loadstep,string name,double (Marker::*func)(unsigned,unsigned)const,unsigned i_,unsigned j_){
		using std::endl;
		fout << "Result \"" << name << "\" \"FEM\" ";
		fout << loadstep << " ";
		fout << "Scalar OnNodes" << endl;
		fout << "Values" <<endl;
		for(unsigned i=0;i<markers_handler->num_markers;i++){
			fout << (i+1) << "\t";
			fout << ((*markers_handler)(i)->*func)(i_,j_) << " ";
			fout << endl;
		}
		fout << "End Values" <<endl;
	};
	template<class Type>
	void gid_result(ostream& fout,unsigned loadstep,string name,Type (Marker::*func)(void)const){
		using std::endl;
		fout << "Result \"" << name << "\" \"FEM\" ";
		fout << loadstep << " ";
		fout << "Scalar OnNodes" << endl;
		fout << "Values" <<endl;
		for(unsigned i=0;i<markers_handler->num_markers;i++){
			fout << (i+1) << "\t";
			fout << ((*markers_handler)(i)->*func)() << " ";
			fout<<endl;
		}
		fout << "End Values" <<endl;
	};

private:
	std::vector<Element_marker*> velements;
	ostringstream marker_file_out;
	Markers_handler*& markers_handler;

	ostringstream m_result;
};
