/*
ver.03.03
file marker.h
*/

#pragma once
#ifdef __INTEL_COMPILER
#include <omp.h>
#endif
#include <vector>
#include <sstream>
#include "element.h"
#include "matrix_algebraic.h"
#include "face.h"
#include "material.h"

//class Markers_handler;

class Marker{
	typedef Elastic<3,MATRIX33> Default_material;
	friend class Markers_handler;
public:
	Marker(void);
	double getWeight(void)const{return weight;};
	double getDisplacement(unsigned i)const{return x[i]-X[i];};
	double getStress(unsigned i,unsigned j)const{
		return sigma(i,j);
	};
	double getStrain(unsigned i,unsigned j)const{
		return F_(i,j);
	};
	double getPotential(void)const{
		return material->potential(&F_);
	};
	int get_parent(void)const{
		return parent_elem;
	};
	void update_constitutive(void){
		material->getC_e(c,&F_);
		material->getSigma(&sigma,&F_);
	};

	Matrix33 F;
	Matrix33 F_;
	Matrix33 F_prev;
	Matrix33 F_prev_;
	Matrix33 sigma;
	double c[3][3][3][3];
	double x[3];
	double x_prev[3];
	double X[3];
private:
	Default_material* material;
	double weight;
	int parent_elem;
	bool isSurface;
};

class Element_marker:virtual public Element{
public:
	Element_marker(){};
	virtual void move_markers(void)=0;
	virtual bool isInside(double x[3])=0;
	virtual bool isInside(Face& face)=0;
	virtual void calcVolume(void)=0;

	virtual void marker_clear(void){
		vmarker.clear();
		vfaces.clear();
	};
	virtual void marker_add(Marker* marker){
		vmarker.push_back(marker);
	};
	void face_add(Face* face){
		vfaces.push_back(face);
	};
	double get_num_marker(void)const{
		return (double)vmarker.size();
	};
	double get_weight(void)const{
		return m_weight;
	};
protected:
	std::vector<Marker*> vmarker;
	std::vector<Face*> vfaces;
	double m_weight;
	double m_marker_Volume;
	double m_marker_volume;
};

class Markers_handler{
	typedef std::ostringstream ostringstream;
	typedef struct{
		Marker* pmarker[3];
		int node_indices[3];
		int elem_index;
	}Connectivity;
public:
	template<class Elem_type,class Material_type>
	Markers_handler(
		ostringstream& filename,
		std::vector<Elem_type*>& velements,
		Material_type* front_material
	){
		for(unsigned i=0;i<velements.size();i++)
			this->velements.push_back(velements[i]);
		read_gid(filename);
		for(unsigned i=0;i<num_markers;i++)
			markers[i].material=front_material;
	};

	~Markers_handler(void);
	void update_element(void);
	unsigned int num_markers;
	Marker* operator()(unsigned i){return &markers[i];};
private:
	void read_gid(ostringstream& filename);
	static const int NODE_IN_ELEM=3;
	static const int DIM=3;
	Marker* markers;
	Face* boundary;
	std::vector<Connectivity*> connectivity;
	std::vector<Element_marker*> velements;
	unsigned int num_faces;
};

#include "elem_gauss.h"
template<const unsigned DIM,const unsigned NUM_NODES,const unsigned NUM_GP,class F_TYPE,const unsigned NUM_FACES>
class Element_gauss_marker:
	virtual public Element_gauss<DIM,NUM_NODES,NUM_GP,F_TYPE>,
	virtual public Element_marker
{
public:
	typedef Element_gauss<DIM,NUM_NODES,NUM_GP,F_TYPE> Parent;
	typedef Elastic<DIM,F_TYPE> Default_material;
	typedef const unsigned FACE_INDEX_ARRAY[NUM_FACES][DIM];

	ELEM_GAUSS_TYPES;

	using Parent::nodes;
	using Parent::material;
	using Parent::gauss_points;
	using Parent::vmatrix_iterator;

	Element_gauss_marker(
		int* node_index,
		Node* raw_nodes,
		Default_material* material,
		GP_WEIGHT_ARRAY&	gp_weight,
		SHAPE_ARRAY&		N,
		DNDXI_ARRAY&		dNdxi,
		GP_COOR_ARRAY&		gp_coor,
		double (*dNdxi_xi)(unsigned,unsigned,double[DIM]),
		double (*N_xi)(unsigned,double[DIM]),
		FACE_INDEX_ARRAY& face_index
	);

	virtual ~Element_gauss_marker(){};

	//from element_gauss
	void update_loadstep(void);
	void update_iteration(void);
	void calcK_c(Sparse& K);
	void calcK_sigma(Sparse& K);
	void calcT(MATRIX& T);

	//from element_marker
	void marker_clear(void);
	void marker_add(Marker* marker);
	void move_markers(void);
	//fixme: general dimension
	bool isInside(double x[3]);
	bool isInside(Face& face);
	void calcVolume(void);

protected:
	void update_markers(void);
	void calcT_unique(MATRIX& T);
	void get_xi(double xi[DIM],double x[DIM]);

	double (*dNdxi_xi)(unsigned,unsigned,double[DIM]);
	double (*N_xi)(unsigned,double[DIM]);

	std::vector<Default_gauss_point*> closest_gp;
	FACE_INDEX_ARRAY &face_index;
	Face m_Face[NUM_FACES];
	double m_weight_all;
};

