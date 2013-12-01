/*
ver.03.03
file tetra1_marker.h
*/

#pragma once
#include "common.h"
#include "tetra1.h"
#include "tetra1_up.hpp"
#include "marker.h"
#include "face.h"

class Tetra1_marker:virtual public Tetra1,virtual public Element_marker{
public:
	Tetra1_marker(
		int* node_index,
		Node* raw_nodes,
		Default_material* material
	);
	virtual ~Tetra1_marker(){};
	//from Tetra1
	void update_loadstep(void);
	void update_iteration(void);

	void calcK(Sparse& K){
		Tetra1::calcK(K);
		calcK_c(K);
		calcK_sigma(K);
#ifdef ALPHA
		calcK_alpha(K);
#endif
#ifdef MASS
		calcK_mass(K);
#endif
	};
	void calcT(MATRIX& T){
		Tetra1::calcT(T);
		calcT_unique(T);
	};
	void calcT_unique(MATRIX& T);

	void calcK_mass(Sparse& K);
	void calcT_mass(MATRIX& T);
	//from Element_marker
	void move_markers(void);
	bool isInside(double x[3]);
	bool isInside(Face& face);
	void calcVolume(void);

	//others
	bool isCut(void){
		if(0<m_weight && m_weight<1.-ABFEM::CONSTANTS::g_eps && vmarker.size())
			return true;
		else
			return false;
	};
	double get_weight(void)const{
		return m_weight;
	};
	double get_volume(void)const{
		return m_weight*m_Volume;
	};

	//issue11
	bool effective_flag;
	unsigned get_effective_flag(void)const{return effective_flag;};

protected:
	void update_markers(void);

private:
	const static double m_geo_eps;
	const static double m_diffuse;
	const static double m_mass_coef;
	const static int m_VP_div=ABFEM::CONSTANTS::g_VP_div;
	const static int m_num_VP=(m_VP_div-3)*(m_VP_div-2)*(m_VP_div-1)/6;
	void calcK_alpha(Sparse& K);
	void calcK_c(Sparse& K);
	void calcK_sigma(Sparse& K);
	void getVolumeCoor(double dest[m_num_nodes],double x[3]);

	const static int m_surface_node_index[m_num_nodes][m_dim];
	Face m_Face[m_num_nodes];
	Node* m_surface_nodes[m_num_nodes][m_dim];
	double m_shape_weight[m_num_nodes];
};

