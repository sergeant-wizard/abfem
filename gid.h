/*
ver.03.03
file gid.h
*/

#pragma once
#include <iterator>
#include <string>
#include <fstream>
#include "prepost.h"
#include "element.h"
#include "node.h"
#include "tetra8.h"

class Node_gid:public Node{
public:
	Node_gid(void){};
	Node_gid(const unsigned dim,unsigned& cnt,std::istream_iterator<double>& it):Node(dim){
		dof=new Dof[dim];
		it++;
		for(int i=0;i<dim;i++){
			dof[i].x=*it++;
			dof[i].rawindex=cnt++;
			dof[i].X0=dof[i].X=dof[i].x;
		}
	}
	~Node_gid(){
		if(dof)
			delete[] dof;
		dof=NULL;
	};
};

class Gid: virtual public PrePostProcessor{
public:
	typedef std::string string;
	typedef std::ofstream ofstream;
	typedef std::ostream ostream;
	typedef std::ostringstream ostringstream;

	template<class Elem_type,class Material_type>
	Gid(
		std::vector<Dof*>& vdof,
		unsigned& num_dof_dirichlet,
		unsigned& num_dof_nond,
		MATRIX& T,
		std::string projectname,
		std::vector<Elem_type*>& elem_buffer,
		Material_type* material,
		const unsigned DIM,
		const unsigned NODE_IN_ELEM,
		const unsigned NUM_GP
	);

	Gid(
		std::vector<Dof*>& vdof,
		unsigned& num_dof_dirichlet,
		unsigned& num_dof_nond,
		MATRIX& T,
		std::string projectname,
		std::vector<Tetra8*>& elem_buffer,
		Tetra8::Default_material* material,
		const unsigned DIM,
		const unsigned NODE_IN_ELEM,
		const unsigned NUM_GP
	);

	virtual ~Gid();
protected:
	template <class Elem_type,class Material_type>
	bool read_gid_file(
		std::vector<Elem_type*>& elem_buffer,
		Material_type* material
	);

	virtual void save_result(unsigned loadstep);
	virtual void post_process(unsigned loadstep);
	void gid_result(
		ostream& fout,unsigned loadstep,string name,
		double (Element_output::*func)(unsigned,unsigned,unsigned)const,unsigned i_,unsigned j_
	);
	template<class OutType>
	void gid_result(
		ostream& fout,unsigned loadstep,string name,OutType (Element_output::*func)(unsigned)const
	);
	template<class OutType>
	void gid_result(
		ostream& fout,unsigned loadstep,string name,OutType (Element_output::*func)(void)const
	);
	template <class type>
	void gid_result_(
		ostream& fout,unsigned loadstep,string name,double (type::*member)(void)const
	);
	template <class type,class outtype>
	void gid_result(
		ostream& fout,unsigned loadstep,string name,
		outtype (type::*member)(void)const,std::vector<type*>& velements
	);
	template <class type>
	void gid_result(
		ostream& fout,unsigned loadstep,string name,type& values,int offset
	);

	ostringstream m_result;

	void get_type(void);
	string gid_elemtype;
	static const string PATH;
private:
	template<class T>
	void safe_delete(T arg){
		if(arg)
			delete[] arg;
	};

};

