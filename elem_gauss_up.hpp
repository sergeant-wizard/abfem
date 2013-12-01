/*
ver.03.03
file elem_gauss_up.hpp
*/

#pragma once
#include "element.h"
#include "elem_gauss.h"

template<const unsigned DIM,const unsigned NUM_NODES,const unsigned NUM_GP,class F_TYPE>
class Element_gauss_UP:
	virtual public Element_gauss<DIM,NUM_NODES,NUM_GP,F_TYPE>,
	virtual public Element_UP
{
protected:
	typedef Element_gauss<DIM,NUM_NODES,NUM_GP,F_TYPE> Parent;
	typedef Incompressible<DIM,F_TYPE> Default_material;
	ELEM_GAUSS_TYPES;

	using Parent::nodes;
	using Parent::material;
	using Parent::gauss_points;

	Element_gauss_UP(
		int* node_index,
		Node* raw_nodes,
		Default_material* material,
		GP_WEIGHT_ARRAY&	gp_weight,
		SHAPE_ARRAY&		N,
		DNDXI_ARRAY&		dNdxi,
		GP_COOR_ARRAY&		m_gp_coor,
		const unsigned num_dof_p
	):
		Parent(
			node_index,raw_nodes,material,gp_weight,N,dNdxi,m_gp_coor
		),
		Element_UP(num_dof_p),
		num_dof_p(num_dof_p)
	{
	};

	virtual ~Element_gauss_UP(){};

	void initialize_sparse(Sparse& K){
		Parent::initialize_sparse(K);
		for(int a=0;a<NUM_NODES;a++)
			for(int i=0;i<DIM;i++){
				if(nodes[a]->dof[i].isDirichlet)
					continue;
				for(unsigned alpha=0;alpha<num_dof_p;alpha++)
					K.add(nodes[a]->dof[i].mi,m_pressure_id+alpha);
			}
	};

	void update_iteration(void){
		Parent::update_F();
		Parent::update_dNdx();
		for(unsigned gp=0;gp<NUM_GP;gp++){
			material->getC_e(m_c_hat[gp],&gauss_points[gp]->F_);
			material->getSigma(&m_sigma_hat[gp],&gauss_points[gp]->F_);
		}
		update_iteration_unique();
	};

	void calcK_UP(Sparse& K){
		for(unsigned alpha=0;alpha<num_dof_p;alpha++)
		for(unsigned a=0;a<NUM_NODES;a++)
		for(unsigned i=0;i<DIM;i++){
			if(nodes[a]->dof[i].isDirichlet)
				continue;
			int row=nodes[a]->dof[i].mi;
			//issue2
			//K.add(row,m_pressure_id,m_dNdx[a][i]*m_Volume);
			double sum=0;
			for(unsigned gpi=0;gpi<NUM_GP;gpi++){
				Default_gauss_point* gp=gauss_points[gpi];
				sum+=
					gp->coor[alpha]*(m_J_rel_[gpi]+1.)*gp->dNdx[a][i]*
					gp->weight*gp->jacob;
			}
			K.add(row,m_pressure_id+alpha,sum);
		}
	};

	void calcT_p(MATRIX& T){
		for(unsigned alpha=0;alpha<num_dof_p;alpha++){
			double sum=0;
			//issue10
			for(unsigned gpi=0;gpi<NUM_GP;gpi++){
				Default_gauss_point* gp=gauss_points[gpi];
				sum+=m_J_rel_[gpi]*gp->coor[alpha]*gp->weight*gp->jacob;
			}
			T(m_pressure_id_+alpha,0)+=sum;
		}
	};

	double get_pressure(unsigned gp)const{
		double sum=0;
		for(unsigned alpha=0;alpha<num_dof_p;alpha++)
			sum+=gauss_points[gp]->coor[alpha]*m_p[alpha];
		return sum;
	};
	double getStress_hat(unsigned alpha,unsigned i,unsigned j)const{
		return m_sigma_hat[alpha](i,j);
	};
private:
	void update_iteration_unique(void){
		for(unsigned gp=0;gp<NUM_GP;gp++){
			m_J_[gp]=
				gauss_points[gp]->F_.I1()+
				gauss_points[gp]->F_.I2()+
				gauss_points[gp]->F_.I3();
			m_J_rel_[gp]=
				gauss_points[gp]->F_rel_.I1()+
				gauss_points[gp]->F_rel_.I2()+
				gauss_points[gp]->F_rel_.I3();
		}

		using ABFEM::MATH::kronecker2;
		for(unsigned gp=0;gp<NUM_GP;gp++){
			double sum=0;
			for(unsigned alpha=0;alpha<num_dof_p;alpha++)
				sum+=m_p[alpha]*gauss_points[gp]->coor[alpha];
			for(int i=0;i<DIM;i++)
			for(int j=0;j<DIM;j++)
			for(int k=0;k<DIM;k++)
			for(int l=0;l<DIM;l++){
				//issue2
				//m_c_p[i][j][k][l]=-m_p*(
					//kronecker2[i][k][j][l]+kronecker2[i][l][j][k]
				//);
				double c_p=sum*(
					kronecker2[i][j][k][l]-(kronecker2[i][k][j][l]+kronecker2[i][l][j][k])
				);
				gauss_points[gp]->c[i][j][k][l]=c_p+m_c_hat[gp][i][j][k][l];
			}
		}
		for(unsigned gp=0;gp<NUM_GP;gp++)
			gauss_points[gp]->sigma.copy(m_sigma_hat[gp]);
		//issue2
		for(unsigned gp=0;gp<NUM_GP;gp++){
			double sum=0;
			for(unsigned alpha=0;alpha<num_dof_p;alpha++)
				sum+=gauss_points[gp]->coor[alpha]*m_p[alpha];
			for(int i=0;i<DIM;i++)
				gauss_points[gp]->sigma(i,i)+=sum*(m_J_rel_[gp]+1.);
		}

	}

	const unsigned num_dof_p;
	double m_J_[NUM_GP];
	double m_J_rel_[NUM_GP];
	double m_c_hat[NUM_GP][DIM][DIM][DIM][DIM];
	Matrix33 m_sigma_hat[NUM_GP];
};

