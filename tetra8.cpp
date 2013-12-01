/*
ver.03.03
file tetra8.cpp
*/

#include "gid.h"

#include "tetra8.h"
#include "gid.hpp"

using std::vector;

Gid::Gid(
	vector<Dof*>& vdof,
	unsigned& num_dof_dirichlet,
	unsigned& num_dof_nond,
	MATRIX& T,
	string projectname,
	vector<Tetra8*>& elem_buffer,
	Tetra8::Default_material* material,
	const unsigned DIM,
	const unsigned NODE_IN_ELEM,
	const unsigned NUM_GP
):PrePostProcessor(vdof,num_dof_dirichlet,num_dof_nond,T,projectname,DIM,NODE_IN_ELEM,NUM_GP)
{
	using std::ostringstream;
	using std::endl;
	read_gid_file<Tetra8,Tetra8::Default_material>(elem_buffer,material);
	for(unsigned i=0;i<elem_buffer.size();i++)
		for(unsigned elem=0;elem<Tetra8::ELEM_IN_ELEM;elem++)
			velements.push_back(elem_buffer[i]->sub(elem));

	ostringstream output_mesh;
	output_mesh << "gid-MAFINS/" <<
		projectname << "-element.gid/" <<
		projectname << "-element.post.msh";
	ofstream fout(output_mesh.str().c_str());
	fout << "Mesh \"post_mesh\" dimension 3 ElemType Tetrahedra Nnode 4" << endl;
	fout << "Coordinates" << endl;
	for(int i=0;i<num_nodes;i++){
		fout << (i+1) << "\t";
		for(int j=0;j<3;j++)
			fout << raw_nodes[i].dof[j].X0 << "\t";
		fout<<endl;
	}
	fout << "end coordinates" << endl;
	fout << "Elements" << endl;
	for(int i=0;i<velements.size();i++){
		fout << (i+1) << "\t";
		for(int j=0;j<4;j++)
			fout << (velements[i]->get_node(j)->dof[0].rawindex/3+1) << "\t";
		fout<<endl;
	}
	fout << "end elements" << endl;

	get_type();
};

const unsigned Tetra8::EDGE[6][2]={
	{0,1},
	{1,2},
	{2,0},
	{0,3},
	{1,3},
	{2,3}
};

const unsigned Tetra8::ELEM_CONNECT[ELEM_IN_ELEM][NODE_SMALL_ELEM]={
	{4,7,6,0},
	{5,8,4,1},
	{5,6,9,2},
	{7,8,9,3},
	{5,8,6,4},
	{8,7,6,4},
	{5,6,8,9},
	{6,7,8,9},
};

const short Tetra8::SUP_TO_SUB[ELEM_IN_ELEM][NODE_IN_ELEM]={
	{ 3,-1,-1,-1, 0,-1, 2, 1,-1,-1},
	{-1, 3,-1,-1, 2, 0,-1,-1, 1,-1},
	{-1,-1, 3,-1,-1, 0, 1,-1,-1, 2},
	{-1,-1,-1, 3,-1,-1,-1, 0, 1, 2},
	{-1,-1,-1,-1, 3, 0, 2,-1, 1,-1},
	{-1,-1,-1,-1, 3,-1, 2, 1, 0,-1},
	{-1,-1,-1,-1,-1, 0, 1,-1, 2, 3},
	{-1,-1,-1,-1,-1,-1, 0, 1, 2, 3}
};

Tetra8::Tetra8(
	int* indices,
	Node* nodes,
	Default_material* material
):
	material(material),
	volume_ratio(1.0),
	Element_UP(1)
{
	for(int i=0;i<NODE_IN_ELEM;i++)
		this->nodes[i]=&nodes[indices[i]];

	for(int elem=0;elem<ELEM_IN_ELEM;elem++){
		int sub_indices[NODE_SMALL_ELEM]={
			indices[ELEM_CONNECT[elem][0]],
			indices[ELEM_CONNECT[elem][1]],
			indices[ELEM_CONNECT[elem][2]],
			indices[ELEM_CONNECT[elem][3]]
		};
		Tetra1_UP* newElem=new Tetra1_UP(
			sub_indices,
			nodes,
			material
		);
		subelements[elem]=newElem;
	}
}

Tetra8::~Tetra8(void){
	for(int i=0;i<ELEM_IN_ELEM;i++){
		delete subelements[i];
	}
}

void Tetra8::initialize_sparse(Sparse& K){
/*for B-bar
	for(int a=0;a<NODE_IN_ELEM;a++)
		for(int b=0;b<NODE_IN_ELEM;b++)
			for(int i=0;i<DIM;i++){
				MATRIXINDEX row=nodes[a]->ci[i].mi;
				for(int j=0;j<DIM;j++){
					MATRIXINDEX col=nodes[b]->ci[j].mi;
					if(col<row || row<0)
						continue;
					global->K.add(row,col);
				}
			}
*/
	for(int elem=0;elem<ELEM_IN_ELEM;elem++)
		subelements[elem]->initialize_sparse(K);
}

//issue12
void Tetra8::calcBodyForce(MATRIX& T){
	for(unsigned elem=0;elem<ELEM_IN_ELEM;elem++)
			subelements[elem]->Tetra1::calcBodyForce(T);
}
		

/*obsolete
double Tetra8::get_J(void){
	double sum=0;
	for(int elem=0;elem<ELEM_IN_ELEM;elem++)
		sum+=log1p(subelements[elem]->get_J_())
			*subelements[elem]->get_Volume();
	return sum;
}

void Tetra8::increment_pressure(double dp){
	m_p+=dp;
	for(int elem=0;elem<ELEM_IN_ELEM;elem++)
		subelements[elem]->get_pressure()+=dp;
}

*/
