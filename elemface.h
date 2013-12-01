/*
ver.03.03
file elemface.h
*/

#pragma once
#include <string>
#include <math.h>
#include "common.h"
#include "tetra1_up.hpp"
#include "sparse.h"

class Elem_face{
public:
	typedef Tetra1_UP INPUT_ELEMENT;
	typedef std::string string;
	typedef std::ifstream ifstream;
	Elem_face(Node* a,Node* b,Node* c,INPUT_ELEMENT* elem,string filename);
	void initialize(Sparse& K);
	void operator+=(Elem_face& right);
	friend bool operator==(Elem_face& left,Elem_face& right);
	virtual void update_iteration(void);
	void calcK(Sparse& K);
	void calcT(MATRIX& T);

protected:
	bool isInterior(void){return elem1;};

	INPUT_ELEMENT* elem0;
	INPUT_ELEMENT* elem1;
	double coef;
	double area;
	double m_delta_h;
	double m_delta_h_weak;
private:
	Node* nodes[3];
	int row;
	int col;
	const bool marker;
};


