/*
ver.03.03
file gid.hpp
*/

#pragma once
#include "gid.h"

//figure out a way not to use these
using std::vector;
using std::endl;

template<class OutType>
void Gid::gid_result(
	ostream& fout,unsigned loadstep,string name,
	OutType (Element_output::*func)(unsigned)const
){
	fout << "Result \"" << name << "\" \"FEM\" ";
	fout << loadstep << " ";
	fout << "Scalar OnGaussPoints gauss_points" << endl;
	fout << "Values" <<endl;
	for(unsigned i=0;i<velements.size();i++){
		fout << (i+1) << "\t";
		for(unsigned gp=0;gp<NUM_GP;gp++)
			fout << (velements[i]->*func)(gp) << endl;
	}
	fout << "End Values" <<endl;
}

template<class OutType>
void Gid::gid_result(
	ostream& fout,unsigned loadstep,string name,
	OutType (Element_output::*func)(void)const
){
	fout << "Result \"" << name << "\" \"FEM\" ";
	fout << loadstep << " ";
	fout << "Scalar OnGaussPoints gauss_points" << endl;
	fout << "Values" <<endl;
	for(unsigned i=0;i<velements.size();i++){
		fout << (i+1) << "\t";
		fout << (velements[i]->*func)() << endl;
	}
	fout << "End Values" <<endl;
}

template<class Elem_type,class Material_type>
Gid::Gid(
	vector<Dof*>& vdof,
	unsigned& num_dof_dirichlet,
	unsigned& num_dof_nond,
	MATRIX& T,
	string projectname,
	vector<Elem_type*>& elem_buffer,
	Material_type* material,
	const unsigned DIM,
	const unsigned NODE_IN_ELEM,
	const unsigned NUM_GP
):PrePostProcessor(vdof,num_dof_dirichlet,num_dof_nond,T,projectname,DIM,NODE_IN_ELEM,NUM_GP)
{
	read_gid_file<Elem_type,Material_type>(elem_buffer,material);
	for(unsigned i=0;i<elem_buffer.size();i++)
		velements.push_back(elem_buffer[i]);

	get_type();
};

template <class Elem_type,class Material_type>
bool Gid::read_gid_file(
	vector<Elem_type*>& elem_buffer,
	Material_type* material
){
	using namespace std;
	ostringstream filename;
	filename << PATH
		<< projectname << "-element.gid/" 
		<< projectname << "-element.dat";

	static const int MAX_CHAR=0x100;
	unsigned raw_index_cnt=0;
	char dump[MAX_CHAR];
	char title[MAX_CHAR];
	int num_elements;


	ifstream fin(filename.str().c_str());
	try{
		if(fin.fail())
			throw filename.str();
	}
	catch(string str){
		cout<< "bad filename: " << str << endl;
		return false;
	}

	fin.getline(dump,MAX_CHAR);
	//title
	fin >> dump >> title;
	for(int i=0;i<3;i++)
		fin.getline(dump,MAX_CHAR);

	//nodes
	fin >> dump >> num_nodes;

	for(int i=0;i<2;i++)
		fin.getline(dump,MAX_CHAR);
	raw_nodes=new Node_gid[num_nodes];
	istream_iterator<double> it(fin);
	for(unsigned i=0;i<num_nodes;i++)
		new(raw_nodes+i) Node_gid(DIM,raw_index_cnt,it);
	fin.clear();
	fin.getline(dump,MAX_CHAR);

	//elements
	fin >> dump >> num_elements;
	for(int i=0;i<2;i++)
		fin.getline(dump,MAX_CHAR);

	istream_iterator<int> it_(fin);
	for(int i=0;i<num_elements;i++){
		int node_index[NODE_IN_ELEM];
		it_++;
		for(int i=0;i<NODE_IN_ELEM;i++)
			node_index[i]=*it_++-1;
		it_++;
		Elem_type *newElem=new Elem_type(
			node_index,raw_nodes,material
		);
		elem_buffer.push_back(newElem);
	}

	//dirichlet condition
	fin.clear();
	for(int i=0;i<2;i++)
		fin.getline(dump,MAX_CHAR);
	num_dof_dirichlet=0;
	while(1){
		int isDirichlet[DIM];
		double dirichlet[DIM];
		int node_index;
		fin>>node_index;
		for(int i=0;i<DIM;i++)
			fin>>isDirichlet[i];
		for(int i=0;i<DIM;i++)
			fin>>dirichlet[i];

		if(fin.fail()){
			fin.clear();
			fin.getline(dump,MAX_CHAR);
			break;
		}

		for(int j=0;j<DIM;j++){
			if(isDirichlet[j]){
				num_dof_dirichlet++;
				raw_nodes[node_index-1].dof[j].isDirichlet=true;
				raw_nodes[node_index-1].dof[j].dirichlet=dirichlet[j];
			}
		}
	}
		
	/*EZ dirichlet
	num_dof_dirichlet=0;
	for(unsigned a=0;a<num_nodes;a++){
		if(fabs(raw_nodes[a].dof[2].X - 1.0) < 1.0E-8){
			for(unsigned i=0;i<3;i++)
				raw_nodes[a].dof[i].isDirichlet=true;
			raw_nodes[a].dof[0].dirichlet=0.0;
			raw_nodes[a].dof[1].dirichlet=0.0;
			raw_nodes[a].dof[2].dirichlet=-0.5;
			num_dof_dirichlet+=3;
		}else if(fabs(raw_nodes[a].dof[2].X - 0.0) < 1.0E-8){
			for(unsigned i=0;i<3;i++){
				raw_nodes[a].dof[i].isDirichlet=true;
				raw_nodes[a].dof[i].dirichlet=0.0;
			}
			num_dof_dirichlet+=3;
		}
	}
	*/

	//sort
	for(int i=0;i<num_nodes;i++){
		for(int j=0;j<DIM;j++){
			if(raw_nodes[i].dof[j].isDirichlet)
				vdof.insert(vdof.begin(),&raw_nodes[i].dof[j]);
			else
				vdof.push_back(&raw_nodes[i].dof[j]);
		}
	}

	raw_to_mi=new int[vdof.size()];
	for(unsigned i=0;i<vdof.size();i++){
		vdof[i]->mi=i-num_dof_dirichlet;
		vdof[i]->fmi=i;
		raw_to_mi[vdof[i]->rawindex]=i;
	}

	num_dof_nond=vdof.size()-num_dof_dirichlet;
	return true;
}

template <class type>
void Gid::gid_result_(
	ostream& fout,unsigned loadstep,string name,
	double (type::*member)(void)const
){
	fout << "Result \"" << name << "\" \"FEM\" ";
	fout << loadstep << " ";
	fout << "Vector OnNodes" << endl;
	fout << "Values" <<endl;
	for(unsigned i=0;i<num_nodes;i++){
		fout << (i+1) << "\t";
		for(unsigned j=0;j<DIM;j++){
			int mi=raw_to_mi[i*DIM+j];
			fout << (vdof[mi]->*member)() << " ";
		}
		fout<<endl;
	}
	fout << "End Values" <<endl;
}

template <class type,class outtype>
void Gid::gid_result(
	ostream& fout,unsigned loadstep,string name,
	outtype (type::*func)(void)const,vector<type*>& velements
){
	fout << "Result \"" << name << "\" \"FEM\" ";
	fout << loadstep << " ";
	fout << "Scalar OnGaussPoints gauss_points" << endl;
	fout << "Values" <<endl;
	for(unsigned i=0;i<velements.size();i++){
		fout << (i+1) << "\t";
		fout << (velements[i]->*func)() << endl;
	}
	fout << "End Values" <<endl;
}

template <class type>
void Gid::gid_result(
	ostream& fout,unsigned loadstep,string name,type& values,int offset=0
){
	fout << "Result \"" << name << "\" \"FEM\" ";
	fout << loadstep << " ";
	fout << "Vector OnNodes" << endl;
	fout << "Values" <<endl;
	for(unsigned int i=0;i<num_nodes;i++){
		fout << (i+1) << "\t";
		for(int j=0;j<DIM;j++){
			int mi=raw_to_mi[i*DIM+j]+offset;
			fout << values[mi] << " ";
		}
		fout<<endl;
	}
	fout << "End Values" <<endl;
}
