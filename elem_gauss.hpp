/*
ver.03.03
file elem_gauss.hpp
*/

#pragma once
#include "elem_gauss.h"

#define TEMPLATE \
template<const unsigned DIM,const unsigned NUM_NODES,const unsigned NUM_GP,class F_TYPE>
#define EG Element_gauss<DIM,NUM_NODES,NUM_GP,F_TYPE>

TEMPLATE
EG::Element_gauss(
	int* node_index,
	Node* raw_nodes,
	Default_material* material,
	GP_WEIGHT_ARRAY&	gp_weight,
	SHAPE_ARRAY&		N,
	DNDXI_ARRAY&		dNdxi,
	GP_COOR_ARRAY&		gp_coor
):
	material(material)
{
	for(unsigned gp=0;gp<NUM_GP;gp++){
		gauss_points[gp]=new Default_gauss_point(
			gp_weight[gp],N[gp],dNdxi[gp],gp_coor[gp]
		);
		gauss_points[gp]->F_    .zero_out();
		gauss_points[gp]->F_rel_.zero_out();
		gauss_points[gp]->F_rel .identity();
	}
	for(unsigned a=0;a<NUM_NODES;a++)
		nodes[a]=&raw_nodes[node_index[a]];

};

TEMPLATE
EG::~Element_gauss(){
	for(unsigned gp=0;gp>NUM_GP;gp++)
		delete gauss_points[gp];
}

TEMPLATE
void EG::initialize(void){
	update_loadstep();
}

TEMPLATE
void EG::initialize_sparse(Sparse& K){
	for(int a=0;a<NUM_NODES;a++)
	for(int b=0;b<NUM_NODES;b++)
		for(int i=0;i<DIM;i++){
			int row=nodes[a]->dof[i].mi;
			if(nodes[a]->dof[i].isDirichlet)
				continue;
			for(int j=0;j<DIM;j++){
				int col=nodes[b]->dof[j].mi;
				if(col<row)
					continue;
				K.add(row,col);
			}
		}
}

TEMPLATE
void EG::update_loadstep(void){
	Matrix33 dxidX[NUM_GP];
	for(unsigned gpi=0;gpi<NUM_GP;gpi++){
		Default_gauss_point* gp=gauss_points[gpi];
		double dXdxi[DIM][DIM];
		for(unsigned i=0;i<DIM;i++)
		for(unsigned j=0;j<DIM;j++){
			double sum=0;
			for(unsigned node=0;node<NUM_NODES;node++)
				sum+=nodes[node]->dof[i].X*gp->dNdxi[node][j];
			dXdxi[i][j]=sum;
		}
		gp->jacob_prev=
			dXdxi[0][0]*dXdxi[1][1]*dXdxi[2][2] +
			dXdxi[0][1]*dXdxi[1][2]*dXdxi[2][0] +
			dXdxi[0][2]*dXdxi[1][0]*dXdxi[2][1] -
			dXdxi[0][2]*dXdxi[1][1]*dXdxi[2][0] -
			dXdxi[0][0]*dXdxi[1][2]*dXdxi[2][1] -
			dXdxi[0][1]*dXdxi[1][0]*dXdxi[2][2];
		double det=1./gp->jacob_prev;

		dxidX[gpi](0,0)=(dXdxi[1][1]*dXdxi[2][2] - dXdxi[1][2]*dXdxi[2][1])*det;
		dxidX[gpi](0,1)=(dXdxi[0][2]*dXdxi[2][1] - dXdxi[0][1]*dXdxi[2][2])*det;
		dxidX[gpi](0,2)=(dXdxi[0][1]*dXdxi[1][2] - dXdxi[0][2]*dXdxi[1][1])*det;
		dxidX[gpi](1,0)=(dXdxi[1][2]*dXdxi[2][0] - dXdxi[1][0]*dXdxi[2][2])*det;
		dxidX[gpi](1,1)=(dXdxi[0][0]*dXdxi[2][2] - dXdxi[0][2]*dXdxi[2][0])*det;
		dxidX[gpi](1,2)=(dXdxi[0][2]*dXdxi[1][0] - dXdxi[0][0]*dXdxi[1][2])*det;
		dxidX[gpi](2,0)=(dXdxi[1][0]*dXdxi[2][1] - dXdxi[1][1]*dXdxi[2][0])*det;
		dxidX[gpi](2,1)=(dXdxi[0][1]*dXdxi[2][0] - dXdxi[0][0]*dXdxi[2][1])*det;
		dxidX[gpi](2,2)=(dXdxi[0][0]*dXdxi[1][1] - dXdxi[0][1]*dXdxi[1][0])*det;
	}

	for(unsigned gp=0;gp<NUM_GP;gp++)
	for(unsigned node=0;node<NUM_NODES;node++)
	for(unsigned i=0;i<DIM;i++){
		double sum=0;
		for(unsigned j=0;j<DIM;j++)
			sum+=gauss_points[gp]->dNdxi[node][j]*dxidX[gp](j,i);

		gauss_points[gp]->dNdX[node][i]=sum;
	}

	for(unsigned gp=0;gp<NUM_GP;gp++){
		gauss_points[gp]->F_prev_.copy(gauss_points[gp]->F_);
		gauss_points[gp]->sigma_prev.copy(gauss_points[gp]->sigma);

		for(unsigned i=0;i<DIM;i++){
			double sum=0;
			for(unsigned a=0;a<NUM_NODES;a++)
				sum+=gauss_points[gp]->N[a]*nodes[a]->dof[i].X;
			gauss_points[gp]->X[i]=sum;
		}
	}
}

TEMPLATE
void EG::update_iteration(void){
	update_F();
	update_Constitutive();
	update_dNdx();
}

TEMPLATE
void EG::calcT(MATRIX& T){
	for(int a=0;a<NUM_NODES;a++)
	for(int i=0;i<DIM;i++){
		double sum=0;
		for(unsigned gpi=0;gpi<NUM_GP;gpi++){
			Default_gauss_point* gp=gauss_points[gpi];
			double sum_=0;
			for(int j=0;j<DIM;j++)
				sum_+=gp->sigma(i,j)*gp->dNdx[a][j];
			sum+=sum_*gp->weight*gp->jacob;
		}
		T(nodes[a]->dof[i].fmi,0)+=sum;
	}

	//issue12
//	for(int a=0;a<NUM_NODES;a++){
//		double sum=0;
//		for(unsigned gp=0;gp<NUM_GP;gp++){
//			sum+=m_N[gp][a]*m_gp_weight[gp]*m_jacob[gp]*0.001;
//		}
//		T(nodes[a]->dof[2].fmi,0)+=sum;
//	}
}

TEMPLATE
void EG::sync_Kmatrix(Sparse& K,bool symmetric){
	vmatrix_iterator.clear();
	for(int a=0;a<NUM_NODES;a++)
	for(int b=0;b<NUM_NODES;b++)
	for(int i=0;i<DIM;i++)
	for(int j=0;j<DIM;j++){
		int row=nodes[a]->dof[i].mi;
		int col=nodes[b]->dof[j].mi;
		if( (symmetric && col<row) || nodes[a]->dof[i].isDirichlet )
			continue;
		vmatrix_iterator.push_back(
			Matrix_iterator(K(row,col),a,b,i,j)
		);
	}
}

//protected
TEMPLATE
void EG::calcK_c(Sparse& K){
	using std::vector;
	for(
		vector<Matrix_iterator>::iterator it=vmatrix_iterator.begin();
		it!=vmatrix_iterator.end();
		it++
	){
		double sum=0.0;
		for(unsigned gpi=0;gpi<NUM_GP;gpi++){
			Default_gauss_point* gp=gauss_points[gpi];
			double sum_=0;
			for(int k=0;k<DIM;k++)
			for(int l=0;l<DIM;l++)
				sum_+=gp->dNdx[it->a][k]*gp->dNdx[it->b][l]*gp->c[it->i][k][it->j][l];
				//sum_+=m_dNdX[gp][it->a][k]*m_dNdX[gp][it->b][l]*m_c[gp][it->i][k][it->j][l];
			sum+=sum_*gp->weight*gp->jacob;
		}
		it->ref+=sum;
	}
}

TEMPLATE
void EG::calcK_sigma(Sparse& K){
	using std::vector;
	for(
		vector<Matrix_iterator>::iterator it=vmatrix_iterator.begin();
		it!=vmatrix_iterator.end();
		it++
	){
		if(it->i != it->j)
			continue;
		double sum=0;
		for(unsigned gpi=0;gpi<NUM_GP;gpi++){
			Default_gauss_point* gp=gauss_points[gpi];
			double sum_=0;
			for(int k=0;k<DIM;k++)
			for(int l=0;l<DIM;l++)
				sum_+=
				gp->sigma(k,l)*
				gp->dNdx[it->a][k]*
				gp->dNdx[it->b][l];
			sum+=sum_*gp->weight*gp->jacob;
		}
		it->ref+=sum;
	}
}

TEMPLATE
void EG::update_F(void){
	for(unsigned gpi=0;gpi<NUM_GP;gpi++){
		Default_gauss_point* gp=gauss_points[gpi];
		for(unsigned i=0;i<DIM;i++)
		for(unsigned j=0;j<DIM;j++){
			double sum=0;
			for(int a=0;a<NUM_NODES;a++)
				sum+=nodes[a]->dof[i].x*gp->dNdX[a][j];
			gp->F_rel(i,j)=sum;
		}

		for(unsigned i=0;i<DIM;i++)
		for(unsigned j=0;j<DIM;j++){
			double sum=0;
			for(int a=0;a<NUM_NODES;a++)
				sum+=nodes[a]->dof[i].u*gp->dNdX[a][j];
			gp->F_rel_(i,j)=sum;
		}
		gp->F_=gp->F_rel_*gp->F_prev_+gp->F_rel_+gp->F_prev_;

		//m_J_=m_F_.I3()+m_F_.I2()+m_F_.I1();
		//m_J_rel_=m_F_rel_.I3()+m_F_rel_.I2()+m_F_rel_.I1();
	}

	for(unsigned gp=0;gp<NUM_GP;gp++)
		gauss_points[gp]->jacob=gauss_points[gp]->jacob_prev*gauss_points[gp]->F_rel.I3();
}

TEMPLATE
void EG::update_Constitutive(void){
	for(unsigned gpi=0;gpi<NUM_GP;gpi++){
		Default_gauss_point* gp=gauss_points[gpi];
		material->getC_e(gp->c,&gp->F_);
		material->getSigma(&gp->sigma,&gp->F_);
	}
}

TEMPLATE
void EG::update_dNdx(void){
	for(unsigned gpi=0;gpi<NUM_GP;gpi++){
		Default_gauss_point* gp=gauss_points[gpi];
		Matrix33 Finv=gp->F_rel.inv();
		for(int node=0;node<NUM_NODES;node++)
		for(int i=0;i<DIM;i++){
				gp->dNdx[node][i]=0;
				for(int I=0;I<DIM;I++)
					gp->dNdx[node][i]+=Finv(I,i)*gp->dNdX[node][I];
		}
	}
}

#undef TEMPLATE
#undef EG
