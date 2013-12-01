/*
ver.03.03
file hexa27.h
*/
#pragma once

#include "elem_gauss.hpp"
#include "elem_gauss_up.hpp"

class Hexa27:virtual public Element_gauss<3,27,27,MATRIX33>{
public:
	Hexa27(
		int* node_index,
		Node* raw_nodes,
		Default_material* material
	):
		Element_gauss<3,27,27,MATRIX33>(
			node_index,raw_nodes,material,
			m_gp_weight,m_N,m_dNdxi,gp_coor
		)
	{
	};
	virtual ~Hexa27(){};

	static double m_dNdxi_xi(unsigned node,unsigned i,double xi[3]);
	static double m_N_xi(unsigned node,double xi[3]);

	static const unsigned m_num_gp=27;
	static const unsigned m_num_nodes=27;
	static const unsigned m_dim=3;
	static const double m_dNdxi[m_num_gp][m_num_nodes][m_dim];
	static const double m_gp_weight[m_num_gp];
	static const double m_N[m_num_gp][m_num_nodes];
	static const double gp_coor[m_num_gp][m_dim+1];
};

class Hexa27_UP:virtual public Hexa27,virtual public Element_gauss_UP<3,27,27,MATRIX33>{
public:
	typedef Incompressible<3,MATRIX33> Default_material;
	Hexa27_UP(
		int* node_index,
		Node* raw_nodes,
		Default_material* material
	):
		Hexa27(node_index,raw_nodes,material),
		Element_gauss<3,27,27,MATRIX33>(
			node_index,raw_nodes,material,
			m_gp_weight,m_N,m_dNdxi,gp_coor
		),
		Element_gauss_UP<3,27,27,MATRIX33>(
			node_index,raw_nodes,material,
			m_gp_weight,m_N,m_dNdxi,gp_coor,4
		),
		Element_UP(4)
	{};
	virtual ~Hexa27_UP(){};
};

#include "marker.h"
#include "elem_gauss_marker.hpp"
class Hexa27_marker:
	virtual public Hexa27,
	virtual public Element_gauss_marker<3,27,27,MATRIX33,12>
{
public:
	Hexa27_marker(
		int* node_index,
		Node* raw_nodes,
		Default_material* material
	):
		Hexa27(node_index,raw_nodes,material),
		Element_gauss<3,27,27,MATRIX33>(
			node_index,raw_nodes,material,
			m_gp_weight,m_N,m_dNdxi,gp_coor
		),
		Element_gauss_marker<3,27,27,MATRIX33,12>(
			node_index,raw_nodes,material,
			m_gp_weight,m_N,m_dNdxi,gp_coor,
			m_dNdxi_xi,
			m_N_xi,
			m_face_index
		)
	{};
	virtual ~Hexa27_marker(){};
private:
	static const unsigned m_face_index[12][3];
};
