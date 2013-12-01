/*
ver.03.03
file elemface.cpp
*/

#include "elemface.h"

Elem_face::Elem_face(Node* n0,Node* n1,Node* n2,Tetra1_UP* elem,string filename):
elem0(elem),elem1(NULL),marker(false)
{
	using ABFEM::MATH::SQR;
	nodes[0]=n0;
	nodes[1]=n1;
	nodes[2]=n2;
	ABFEM::SORT(3,nodes);

	ifstream fin(filename.c_str());
	//issue5
	try{
		if(fin.fail()){
			std::cout << "elemface failed to open conf file" << std::endl;
			throw 0;
		}
	}
	catch(int){return;}

	double delta,h,weak_coef;
	fin >> delta;
	//fin >> h;
	h=sqrt(
		SQR(n0->dof[0].X-n1->dof[0].X)+
		SQR(n0->dof[1].X-n1->dof[1].X)+
		SQR(n0->dof[2].X-n1->dof[2].X)
	);
	fin >> weak_coef;
	m_delta_h=delta*h;
	m_delta_h_weak=m_delta_h*weak_coef;

};

void Elem_face::initialize(Sparse& K){
	if(!isInterior())
		return;
	row=elem0->get_pressureID();
	col=elem1->get_pressureID();
	if(col<row){
		int tmp=row;
		row=col;
		col=tmp;
	}
	K.add(row,col);
};

void Elem_face::operator+=(Elem_face& right){
	elem1=right.elem0;
};

bool operator==(Elem_face& left,Elem_face& right){
	for(int i=0;i<3;i++)
		if(left.nodes[i]!=right.nodes[i])
			return false;
	return true;
};

void Elem_face::update_iteration(void){
	using namespace ABFEM::MATH;
	if(!isInterior())
		return;
	double a0=
		nodes[0]->dof[1].x*nodes[1]->dof[2].x-nodes[0]->dof[2].x*nodes[1]->dof[1].x+
		nodes[1]->dof[1].x*nodes[2]->dof[2].x-nodes[1]->dof[2].x*nodes[2]->dof[1].x+
		nodes[2]->dof[1].x*nodes[0]->dof[2].x-nodes[2]->dof[2].x*nodes[0]->dof[1].x;
	double a1=
		nodes[0]->dof[2].x*nodes[1]->dof[0].x-nodes[0]->dof[0].x*nodes[1]->dof[2].x+
		nodes[1]->dof[2].x*nodes[2]->dof[0].x-nodes[1]->dof[0].x*nodes[2]->dof[2].x+
		nodes[2]->dof[2].x*nodes[0]->dof[0].x-nodes[2]->dof[0].x*nodes[0]->dof[2].x;
	double a2=
		nodes[0]->dof[0].x*nodes[1]->dof[1].x-nodes[0]->dof[1].x*nodes[1]->dof[0].x+
		nodes[1]->dof[0].x*nodes[2]->dof[1].x-nodes[1]->dof[1].x*nodes[2]->dof[0].x+
		nodes[2]->dof[0].x*nodes[0]->dof[1].x-nodes[2]->dof[1].x*nodes[0]->dof[0].x;
	area=.5*sqrt(SQR(a0)+SQR(a1)+SQR(a2));

	coef=area*m_delta_h;

	return;
};

void Elem_face::calcK(Sparse& K){
	if(!isInterior())
		return;

	K.add(row,row,coef);
	K.add(col,col,coef);
	K.add(row,col,-coef);
}

void Elem_face::calcT(MATRIX& T){
	if(!isInterior())
		return;
	T(elem0->get_pressureID_(),0)+=coef*(
		elem0->get_pressure(0) - elem1->get_pressure(0)
	);
	T(elem1->get_pressureID_(),0)+=coef*(
		elem1->get_pressure(0) - elem0->get_pressure(0)
	);

}

