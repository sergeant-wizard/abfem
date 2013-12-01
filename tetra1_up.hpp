/*
ver.03.03
file tetra1_up.hpp
*/

#pragma once
#include "tetra1.h"

/*
difference with ln(J) and J-1...
1. sigma_p
2. T_p
3. K_UP
4. c_p
*/

class Tetra1_UP:virtual public Tetra1,virtual public Element_UP{
public:
	typedef Incompressible<3,MATRIX33> Default_material;
	Tetra1_UP(
		int* node_index,
		Node* raw_nodes,
		Default_material* material
	):
		Tetra1(node_index,raw_nodes,material),
		Element_UP(1)
	{};
	~Tetra1_UP(){};

	void initialize_sparse(Sparse& K){
		Tetra1::initialize_sparse(K);
		for(int a=0;a<m_num_nodes;a++)
			for(int i=0;i<m_dim;i++){
				if(nodes[a]->dof[i].isDirichlet)
					continue;
				K.add(nodes[a]->dof[i].mi,m_pressure_id);
			}
	};
	void update_loadstep(void){
		Tetra1::update_loadstep();
		m_p[0]=0;
	};

	virtual void update_iteration(void){
		Tetra1::update_F();
		Tetra1::update_dNdx();
		material->getC_e(m_c_hat,&m_F_);
		material->getSigma(&m_sigma_hat,&m_F_);
		update_iteration_unique();
	};

	void calcK_UP(Sparse& K){
		for(unsigned a=0;a<m_num_nodes;a++)
		for(unsigned i=0;i<m_dim;i++){
			if(nodes[a]->dof[i].isDirichlet)
				continue;
			int row=nodes[a]->dof[i].mi;
			//issue2
			//K.add(row,m_pressure_id,m_dNdx[a][i]*m_volume);
			K.add(row,m_pressure_id,m_dNdx[a][i]*m_volume*(1.+m_J_rel_));
		}
	};
	void calcT_p(MATRIX& T){
		//issue2
		//T(m_pressure_id_,0)+=log1p(m_J_)*m_volume;
		//T(m_pressure_id_,0)+=log1p(m_J_rel_)*m_volume;
		T(m_pressure_id_,0)+=m_J_rel_*m_volume;
	};

	double get_pressure(unsigned)const{
		return m_p[0];
	};
	double getStress_hat(unsigned,unsigned i,unsigned j)const{
		return m_sigma_hat(i,j);
	};
protected:
	void update_iteration_unique(void){
		using ABFEM::MATH::kronecker2;
		for(int i=0;i<m_dim;i++)
		for(int j=0;j<m_dim;j++)
		for(int k=0;k<m_dim;k++)
		for(int l=0;l<m_dim;l++){
			//issue2
			/*
			m_c_p[i][j][k][l]=-m_p*(
				kronecker2[i][k][j][l]+kronecker2[i][l][j][k]
			);
			*/
			//wth we should need J
			m_c_p[i][j][k][l]=m_p[0]*(
				kronecker2[i][j][k][l]-(kronecker2[i][k][j][l]+kronecker2[i][l][j][k])
			);
			m_c[i][j][k][l]=m_c_p[i][j][k][l]+m_c_hat[i][j][k][l];
		}

		m_sigma.copy(m_sigma_hat);
		//issue2
		for(int i=0;i<m_dim;i++)
			//m_sigma(i,i)+=m_p;
			m_sigma(i,i)+=m_p[0]*(1.+m_J_rel_);
	};

	Matrix33 m_sigma_hat;
	double m_c_hat[m_dim][m_dim][m_dim][m_dim];

private:
	double m_c_p[m_dim][m_dim][m_dim][m_dim];
};

