/*
ver.03.03
file node.h
*/

#pragma once

struct Dof{
	double dirichlet;
	unsigned rawindex;
	int mi;
	unsigned fmi;
	bool isDirichlet;
	double X0;
	double X;
	double x;
	double u;
	double get_displacement(void)const{
		return x-X0;
	};
	Dof():
	    isDirichlet(false),
		u(0)
		{};
};

class Node{
public:
	Node(void):dof(NULL),dim(1){};
	Node(const int dim):dim(dim){
		dof=new Dof[dim];
	};
	virtual ~Node(){
		if(dof)
			delete[] dof;
		dof=NULL;
	};
	Dof* dof;
protected:
	const int dim;
};


