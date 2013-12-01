/*
ver.03.03
file elem_gauss_marker.hpp
*/

#pragma once
#include "marker.h"

#define EGM Element_gauss_marker<DIM,NUM_NODES,NUM_GP,F_TYPE,NUM_FACES>
#define TEMPLATE template<const unsigned DIM,const unsigned NUM_NODES,const unsigned NUM_GP,class F_TYPE,const unsigned NUM_FACES>

TEMPLATE
EGM::Element_gauss_marker(
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
):
	Parent(node_index,raw_nodes,material,gp_weight,N,dNdxi,gp_coor),
	dNdxi_xi(dNdxi_xi),
	N_xi(N_xi),
	face_index(face_index)
{
	m_weight_all=0;
	for(unsigned gp=0;gp<NUM_GP;gp++)
		m_weight_all+=gp_weight[gp];
}

TEMPLATE
void EGM::update_loadstep(void){
	Parent::update_loadstep();
	for(unsigned m=0;m<vmarker.size();m++){
		vmarker[m]->F_prev.copy(vmarker[m]->F);
		vmarker[m]->F_prev_.copy(vmarker[m]->F_);
		for(unsigned i=0;i<DIM;i++)
			vmarker[m]->x_prev[i]=vmarker[m]->x[i];
	}

	for(unsigned f=0;f<NUM_FACES;f++)
		for(unsigned i=0;i<DIM;i++)
		for(unsigned j=0;j<DIM;j++)
			m_Face[f].initialize(
				nodes[face_index[f][i]]->dof[j].X
			);
}

TEMPLATE
void EGM::update_iteration(void){
	Parent::update_F();
	Parent::update_dNdx();
	for(unsigned gpi=0;gpi<NUM_GP;gpi++){
		Default_gauss_point* gp=gauss_points[gpi];
		material->getC_e(gp->c,&gp->F_rel_);
		material->getSigma(&gp->sigma,&gp->F_rel_);
	}
	update_markers();
}

TEMPLATE
void EGM::calcK_c(Sparse& K){
	
	using std::vector;
	for(typename vector<typename Parent::Matrix_iterator>::iterator it=vmatrix_iterator.begin();it!=vmatrix_iterator.end();it++){
		if(!vmarker.size() || !m_weight)
			continue;

		double numerator=0;
		double denominator=0;
		for(unsigned m=0;m<vmarker.size();m++){
			Default_gauss_point* gp=closest_gp[m];
			denominator+=vmarker[m]->getWeight();

			double sum=0;
			for(int k=0;k<DIM;k++)
			for(int l=0;l<DIM;l++)
				sum+=
					gp->dNdx[it->a][k]*
					gp->dNdx[it->b][l]*
					vmarker[m]->c[it->i][k][it->j][l];
			numerator+=
				sum*
				vmarker[m]->getWeight()*
				closest_gp[m]->weight*
				closest_gp[m]->jacob;
		}
		it->ref+=numerator/denominator;
	}
}

TEMPLATE
void EGM::calcK_sigma(Sparse& K){
	using std::vector;
	for(typename vector<typename Parent::Matrix_iterator>::iterator it=vmatrix_iterator.begin();it!=vmatrix_iterator.end();it++){
		if(it->i!=it->j)
			continue;
		if(!vmarker.size() || !m_weight)
			continue;

		double numerator=0;
		double denominator=0;
		for(unsigned m=0;m<vmarker.size();m++){
			Default_gauss_point* gp=closest_gp[m];
			denominator+=vmarker[m]->getWeight();

			double sum=0;
			for(int k=0;k<DIM;k++)
			for(int l=0;l<DIM;l++)
				sum+=
					vmarker[m]->sigma(k,l)*
					gp->dNdx[it->a][k]*
					gp->dNdx[it->b][l];
			numerator+=
				sum*
				vmarker[m]->getWeight()*
				closest_gp[m]->weight*
				closest_gp[m]->jacob;
		}
		it->ref+=numerator/denominator;
	}
}

TEMPLATE
void EGM::calcT(MATRIX& T){
	Parent::calcT(T);
	calcT_unique(T);
}

TEMPLATE
void EGM::calcT_unique(MATRIX& T){

	for(int a=0;a<NUM_NODES;a++)
	for(int i=0;i<DIM;i++){
		if(!vmarker.size() || !m_weight){
			continue;
		}

		double numerator=0;
		double denominator=0;
		for(unsigned m=0;m<vmarker.size();m++){
			Default_gauss_point* gp=closest_gp[m];
			denominator+=vmarker[m]->getWeight();
			double sum=0;
			for(int j=0;j<DIM;j++)
				sum+=
					vmarker[m]->sigma(i,j)*
					gp->dNdx[a][j];
			numerator+=
				sum*
				vmarker[m]->getWeight()*
				closest_gp[m]->weight*
				closest_gp[m]->jacob;
		}
		T(nodes[a]->dof[i].fmi,0)+=numerator/denominator;
	}
}

TEMPLATE
void EGM::marker_clear(void){
	Element_marker::marker_clear();
	closest_gp.clear();
}

TEMPLATE
void EGM::marker_add(Marker* marker){
	Element_marker::marker_add(marker);
	for(unsigned m=0;m<vmarker.size();m++)
		closest_gp.push_back(gauss_points[26]);
		//fixme: always center for now
}

TEMPLATE
void EGM::move_markers(void){
	for(unsigned m=0;m<vmarker.size();m++){
		double xi[DIM];
		get_xi(xi,vmarker[m]->x_prev);

		for(unsigned i=0;i<DIM;i++){
			double sum=0;
			for(unsigned a=0;a<NUM_NODES;a++)
				sum+=nodes[a]->dof[i].u*(*N_xi)(a,xi);

			vmarker[m]->x[i]=
				vmarker[m]->x_prev[i]+sum;
		}
	}
}

TEMPLATE
bool EGM::isInside(double x[3]){
	using ABFEM::CONSTANTS::g_geo_eps;
	double xi[3];
	get_xi(xi,x);

	for(unsigned i=0;i<3;i++)
		if(xi[i]<-1.-g_geo_eps || 1.+g_geo_eps<xi[i])
			return false;
	return true;
}

TEMPLATE
bool EGM::isInside(Face& face){
	for(int i=0;i<DIM;i++)
		if(isInside(face.getVertice(i)))
			return true;

	for(int i=0;i<NUM_FACES;i++){
		if(Face::Cross(face,m_Face[i])){
			return true;
		}
	}
	return false;
}

TEMPLATE
void EGM::calcVolume(void){
	if(!vfaces.size()){
		if(vmarker.size()){
			m_weight=m_weight_all;
		}else{
			m_weight=0.0;
		}
		return;
	}

	double weight_sum=0;
	for(int gp=0;gp<NUM_GP;gp++){
		bool flag=true;
		for(unsigned int face=0;face<vfaces.size();face++)
			if(vfaces[face]->distance(gauss_points[gp]->X)<0)
				flag=false;

		if(flag){
			weight_sum+=gauss_points[gp]->weight;
		}
	}

	//issue3
	//if(!weight_sum)
		//weight_sum=1./(double)m_num_VP;

#ifdef MAX_VOLUME
	if(0<weight_sum && weight_sum < 1.-1E-8)
		weight_sum=m_weight_all;
#endif

	//m_marker_Volume=m_Volume*weight_sum;
	m_weight=weight_sum;
}

TEMPLATE
void EGM::update_markers(void){
	for(unsigned m=0;m<vmarker.size();m++){
		Default_gauss_point* gp=closest_gp[m];
		//check: are we really not using marker.F?
		//vmarker[m]->F=m_F_rel[gp]*vmarker[m]->F_prev;
		vmarker[m]->F_=
			gp->F_rel_*vmarker[m]->F_prev_+
			gp->F_rel_+
			vmarker[m]->F_prev_;


		vmarker[m]->update_constitutive();
	};
}

TEMPLATE
void EGM::get_xi(double xi[DIM],double x[DIM]){
	static const double tol=1.0E-8;
	static const unsigned max_iteration=16;
	double xi_n[DIM]={0};
	double f_n[DIM];
	double K_inv[DIM][DIM];
	unsigned iteration_cnt=0;

	do{
		double norm=0;
		for(unsigned I=0;I<DIM;I++){
			double sum=0;
			for(unsigned a=0;a<NUM_NODES;a++)
				sum+=nodes[a]->dof[I].X*(*N_xi)(a,xi_n);

			f_n[I]=x[I]-sum;
			norm+=fabs(f_n[I]);
		}

		if(norm<tol) 
			break;

		if(iteration_cnt++>max_iteration){
			std::cout << "get_xi failed" << std::endl;
			return;
		}
			
		switch(DIM){
			case 2:{//this hasn't been tested yet
				double K[2][2];
				for(unsigned I=0;I<2;I++)
				for(unsigned i=0;i<2;i++){
					double sum=0;
					for(unsigned a=0;a<NUM_NODES;a++)
						sum+=nodes[a]->dof[I].X*(*dNdxi_xi)(a,i,xi_n);
					K[I][i]=sum;
				}
				double det=1./(K[0][0]*K[1][1]-K[0][1]*K[1][0]);
				K_inv[0][0]= K[1][1]*det;
				K_inv[0][1]=-K[0][1]*det;
				K_inv[1][0]=-K[1][0]*det;
				K_inv[1][1]= K[0][0]*det;
				break;
		   }default:{
				Matrix33 K;
				for(unsigned I=0;I<3;I++)
				for(unsigned i=0;i<3;i++){
					double sum=0;
					for(unsigned a=0;a<NUM_NODES;a++)
						sum+=nodes[a]->dof[I].X*(*dNdxi_xi)(a,i,xi_n);
					K(I,i)=sum;
				}

				Matrix33 K_inv33=K.inv();
				for(unsigned i=0;i<3;i++)
				for(unsigned j=0;j<3;j++)
					K_inv[i][j]=K_inv33(i,j);
				break;
			}
		}

		for(unsigned i=0;i<DIM;i++){
			double delta_xi=0;
			for(unsigned I=0;I<DIM;I++)
				delta_xi+=K_inv[i][I]*f_n[I];
			xi_n[i]+=delta_xi;
		}
	}while(1);

	for(unsigned i=0;i<DIM;i++)
		xi[i]=xi_n[i];
}

#undef EGM
#undef TEMPLATE
