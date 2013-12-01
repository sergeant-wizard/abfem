/*
ver.03.03
file elem_gauss.h
*/

#pragma once

#include <vector>
#include "common.h"
#include "element.h"
#include "material.h"
#include "matrix_algebraic.h"
#include "node.h"
#include "sparse.h"


#define ELEM_GAUSS_TYPES \
	typedef typename Parent::Default_matrix Default_matrix; \
	typedef typename Parent::Default_gauss_point Default_gauss_point; \
	typedef typename Parent::GP_WEIGHT_ARRAY GP_WEIGHT_ARRAY; \
	typedef typename Parent::SHAPE_ARRAY SHAPE_ARRAY; \
	typedef typename Parent::DNDXI_ARRAY DNDXI_ARRAY; \
	typedef typename Parent::GP_COOR_ARRAY GP_COOR_ARRAY; \

#pragma once

template<const unsigned DIM,const unsigned NUM_NODES,class Default_matrix,class REAL=double>
class Gauss_point{
public:
	typedef const double GP_WEIGHT_ARRAY;
	typedef const double SHAPE_ARRAY[NUM_NODES];
	typedef const double DNDXI_ARRAY[NUM_NODES][DIM];
	typedef const double GP_COOR_ARRAY[DIM+1];

	Gauss_point(
		GP_WEIGHT_ARRAY&	weight,
		SHAPE_ARRAY&		N,
		DNDXI_ARRAY&		dNdxi,
		GP_COOR_ARRAY&		coor
	):
		weight(weight),N(N),dNdxi(dNdxi),coor(coor)
	{};

	virtual ~Gauss_point(){};

	Default_matrix F_;
	Default_matrix F_rel_;
	Default_matrix F_rel;
	Default_matrix F_prev_;
	Default_matrix sigma;
	Default_matrix sigma_prev;
	REAL jacob;
	REAL jacob_prev;
	REAL dNdX[NUM_NODES][DIM];
	REAL dNdx[NUM_NODES][DIM];
	REAL c[DIM][DIM][DIM][DIM];

	GP_WEIGHT_ARRAY&	weight;
	SHAPE_ARRAY&		N;
	DNDXI_ARRAY&		dNdxi;
	GP_COOR_ARRAY&		coor;

	REAL X[DIM];
};

class Element_gauss_base:virtual public Element,virtual public Element_output{
protected:
	virtual ~Element_gauss_base(){};
	
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
};

template<const unsigned DIM,const unsigned NUM_NODES,const unsigned NUM_GP,class F_TYPE>
class Element_gauss:virtual public Element_gauss_base{
protected:
	typedef Elastic<DIM,MATRIX33> Default_material;

	typedef const double GP_WEIGHT_ARRAY[NUM_GP];
	typedef const double SHAPE_ARRAY[NUM_GP][NUM_NODES];
	typedef const double DNDXI_ARRAY[NUM_GP][NUM_NODES][DIM];
	typedef const double GP_COOR_ARRAY[NUM_GP][DIM+1];
	typedef Matrix33 Default_matrix;
	typedef Gauss_point<DIM,NUM_NODES,Default_matrix> Default_gauss_point;

	Element_gauss(
		int* node_index,
		Node* raw_nodes,
		Default_material* material,
		GP_WEIGHT_ARRAY&	gp_weight,
		SHAPE_ARRAY&		N,
		DNDXI_ARRAY&		dNdxi,
		GP_COOR_ARRAY&		gp_coor
	);
	virtual ~Element_gauss();

	virtual void calcK_c(Sparse& K);
	virtual void calcK_sigma(Sparse& K);
	void update_F(void);
	void update_Constitutive(void);
	void update_dNdx(void);

	Default_material* material;
	Node* nodes[NUM_NODES];
	Default_gauss_point* gauss_points[NUM_GP];
public:
	virtual void initialize(void);
	virtual void initialize_sparse(Sparse& K);
	virtual void update_loadstep(void);
	virtual void update_iteration(void);
	virtual void calcK(Sparse& K){
		calcK_c(K);
		calcK_sigma(K);
	};
	virtual void calcT(MATRIX& T);
	void sync_Kmatrix(Sparse& K,bool symmetric=true);
	
	//from element_output
	double getStress(unsigned gp,unsigned i,unsigned j)const{
		return gauss_points[gp]->sigma(i,j);
	};
	double get_F_(unsigned gp,unsigned i,unsigned j)const{
		return gauss_points[gp]->F_(i,j);
	};
	double get_J_(unsigned gp)const{
		const Default_matrix& F=gauss_points[gp]->F_;
		return F.I1()+F.I2()+F.I3();
	};
	double get_pressure(unsigned gp)const{
		return gauss_points[gp]->sigma.I1()/3.0;
	};
	double get_potential(unsigned gp)const{
		return material->potential(&gauss_points[gp]->F_);
	};
	Node*& get_node(unsigned node){return nodes[node];};
};

