/*
ver.03.03
file tetra1.h
*/

#pragma once
#include "common.h"
#include "element.h"
#include "material.h"
#include "matrix_algebraic.h"

class Tetra1: virtual public Element, virtual public Element_output{
public:
	typedef Elastic<3,MATRIX33> Default_material;
	Tetra1(
		int* node_index,
		Node* raw_nodes,
		Default_material* material
	);
	virtual ~Tetra1();

	//from Element
	virtual void initialize(void);
	virtual void initialize_sparse(Sparse& K);
	virtual void update_loadstep(void);
	virtual void update_iteration(void);
	void calcK(Sparse& K){
		calcK_c(K);
		calcK_sigma(K);
	};
	virtual void calcT(MATRIX& T);
	//issue12
	void calcBodyForce(MATRIX& T);
	Node*& get_node(unsigned node){
		return nodes[node];
	};
	void sync_Kmatrix(Sparse& K,bool symmetric=true);

	//from Element_output
	double getStress(unsigned,unsigned i,unsigned j)const{
		return m_sigma(i,j);
	};
	double get_F_(unsigned,unsigned i,unsigned j)const{
		return m_F_(i,j);
	};
	double get_J_(unsigned)const{
		return m_J_;
	};
	virtual double get_pressure(unsigned)const{
		return m_sigma.I1()/3.0;
	};
	virtual double get_potential(unsigned)const{
		return material->potential(&m_F_);
	};
	virtual double get_volume(void)const{
		return m_Volume;
	};

	//others
	double get_dndx(unsigned node,unsigned comp)const{
		return m_dNdx[node][comp];
	};
	double get_Volume(void){
		return m_Volume;
	};
protected:
	void calcT_p(MATRIX& T);
	void update_F(void);
	void update_Constitutive(void);
	void update_dNdx(void);

	const static unsigned m_num_nodes=4;
	const static unsigned m_dim=3;

	Default_material* material;

	Matrix33 m_sigma;
	Matrix33 m_sigma_prev;

	Matrix33 m_F;
	Matrix33 m_F_prev;
	Matrix33 m_F_rel;
	
	Matrix33 m_F_;
	Matrix33 m_F_prev_;
	Matrix33 m_F_rel_;

	double m_J;
	double m_J_;
	double m_J_rel_;
	double m_Volume;
	double m_volume;
	double m_dNdX[m_num_nodes][m_dim];
	double m_dNdx[m_num_nodes][m_dim];

	double m_c[m_dim][m_dim][m_dim][m_dim];

	Node* nodes[m_num_nodes];

	//iterators
	struct Matrix_iterator{
		Matrix_iterator(double& ref_,unsigned a_,unsigned b_,unsigned i_,unsigned j_):
		ref(ref_),a(a_),b(b_),i(i_),j(j_){};
		void operator=(Matrix_iterator& right){
			ref=right.ref;
			a=right.a;
			b=right.b;
			i=right.i;
			j=right.j;
		}
		double& ref;
		unsigned a,b,i,j;
	};
	std::vector<Matrix_iterator> vmatrix_iterator;

private:
	void calcK_c(Sparse& K);
	void calcK_sigma(Sparse& K);

	int m_isNeumannSurface[m_num_nodes];

	void debug_(void);
};

