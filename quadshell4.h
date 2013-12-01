#pragma once

#include "elem_gauss.hpp"

class QuadShell4: virtual public Element_gauss<3,4,4,MATRIX33>{
public:
	QuadShell4(
		int* node_index,
		Node* raw_nodes,
		Default_material* material
	):
		Element_gauss<3,4,4,MATRIX33>(
			node_index,raw_nodes,material,
			m_gp_weight,m_N,m_dNdxi,gp_coor
		)
	{
	};
	virtual ~QuadShell4(){};

	static const unsigned m_num_dof=3;
	double Bmat[m_num_dof]

	static const unsigned m_dim=3;
	static const unsigned m_num_gp=4;
	static const unsigned m_num_nodes=4;

	static const double m_gp_weight[m_num_gp];
	static const double m_N[m_num_gp][m_num_nodes];
	static const double m_dNdxi[m_num_gp][m_num_nodes][m_dim];
	static const double gp_coor[m_num_gp][m_dim+1];
};



