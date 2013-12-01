/*
ver.03.03
file tetra1_marker.cpp
*/

#include "tetra1_marker.h"

const double Tetra1_marker::m_geo_eps=ABFEM::CONSTANTS::g_geo_eps;
const double Tetra1_marker::m_diffuse=ABFEM::CONSTANTS::g_diffuse;
const double Tetra1_marker::m_mass_coef=ABFEM::CONSTANTS::g_mass_coef;

Tetra1_marker::Tetra1_marker(
	int* node_index,
	Node* raw_nodes,
	Default_material* material
):
	Tetra1(node_index,raw_nodes,material),
//issue11
	effective_flag(false)
{
	static const int surface_node_index[4][3]={
		{1,2,3},
		{0,3,2},
		{3,0,1},
		{2,1,0}
	};
	for(unsigned i=0;i<m_num_nodes;i++)
	for(unsigned j=0;j<m_dim;j++)
		m_surface_nodes[i][j]=nodes[surface_node_index[i][j]];

};

void Tetra1_marker::update_loadstep(void){
	Tetra1::update_loadstep();
	for(unsigned i=0;i<vmarker.size();i++){
		vmarker[i]->F_prev.copy(vmarker[i]->F);
		vmarker[i]->F_prev_.copy(vmarker[i]->F_);
	}
	for(unsigned i=0;i<m_num_nodes;i++)
	for(unsigned j=0;j<m_dim;j++)
	for(unsigned k=0;k<m_dim;k++)
		m_Face[i].initialize(
			m_surface_nodes[i][j]->dof[k].X
		);
};

void Tetra1_marker::update_iteration(void){
	Tetra1::update_F();
	Tetra1::update_dNdx();
	material->getC_e(m_c,&m_F_rel_);
	material->getSigma(&m_sigma,&m_F_rel_);
	update_markers();
	m_marker_volume=m_marker_Volume*m_volume/m_Volume;
};

void Tetra1_marker::update_markers(void){
	for(unsigned int m=0;m<vmarker.size();m++){
		vmarker[m]->F=m_F_rel*vmarker[m]->F_prev;
		vmarker[m]->F_=m_F_rel_*vmarker[m]->F_prev_+m_F_rel_+vmarker[m]->F_prev_;
		vmarker[m]->update_constitutive();
	}
}

void Tetra1_marker::move_markers(void){
	for(unsigned int i=0;i<vmarker.size();i++){
		double vol_coor[m_num_nodes]={0};
		getVolumeCoor(vol_coor,vmarker[i]->x);
		for(int j=0;j<m_dim;j++){
			double sum=0;
			for(int node=0;node<m_num_nodes;node++)
				sum+=nodes[node]->dof[j].x*vol_coor[node];
			vmarker[i]->x[j]=sum;
		}
	}
}

bool Tetra1_marker::isInside(double x[3]){
	double vol_coor[4]={0};
	getVolumeCoor(vol_coor,x);
	for(int i=0;i<4;i++)
		//strict
		//if(vol_coor[i]<0.+m_geo_eps || 1.-m_geo_eps<vol_coor[i])
		//loose
		if(vol_coor[i]<0.-m_geo_eps || 1.+m_geo_eps<vol_coor[i])
		//exact
		//if(vol_coor[i]<0.|| 1.<vol_coor[i])
			return false;
	return true;
}

bool Tetra1_marker::isInside(Face& face){
	using ABFEM::MATH::SQR;
#ifdef FAST_ISINSIDE
	double* face0=face.getVertice(0);
	if(
		SQR(face0[0]-nodes[0]->dof[0].X)+
		SQR(face0[1]-nodes[0]->dof[1].X)+
		SQR(face0[2]-nodes[0]->dof[2].X) > SQR(0.5)
	){
		return false;
	}
#endif

	for(int i=0;i<3;i++)
		if(isInside(face.getVertice(i)))
			return true;

	for(int i=0;i<4;i++){
		if(Face::Cross(face,m_Face[i])){
			return true;
		}
	}

	return false;
}

void Tetra1_marker::calcVolume(void){
//VP12
/*
	static const int num_VP=12;
	static const struct{
		double coef[4];
		double weight;
	}volume_point[num_VP]={
		{.625 ,.125 ,.125 ,.125 ,.125 },
		{.125 ,.625 ,.125 ,.125 ,.125 },
		{.125 ,.125 ,.625 ,.125 ,.125 },
		{.125 ,.125 ,.125 ,.625 ,.125 },

		{.4375,.1875,.1875,.1875,.0625},
		{.1875,.4375,.1875,.1875,.0625},
		{.1875,.1875,.4375,.1875,.0625},
		{.1875,.1875,.1875,.4375,.0625},

		{.0625,.3125,.3125,.3125,.0625},
		{.3125,.0625,.3125,.3125,.0625},
		{.3125,.3125,.0625,.3125,.0625},
		{.3125,.3125,.3125,.0625,.0625}
	};
*/
//VP32
/*
	static const int num_VP=32;
	static const struct{
		double coef[4];
		double weight;
	}volume_point[num_VP]={
		{9./12.,1./12.,1./12.,1./12.,1./27.},
		{1./12.,9./12.,1./12.,1./12.,1./27.},
		{1./12.,1./12.,9./12.,1./12.,1./27.},
		{1./12.,1./12.,1./12.,9./12.,1./27.},

		{27./48., 7./48., 7./48., 7./48.,5./108.},
		{ 7./48.,27./48., 7./48., 7./48.,5./108.},
		{ 7./48., 7./48.,27./48., 7./48.,5./108.},
		{ 7./48., 7./48., 7./48.,27./48.,5./108.},

		{23./48., 3./48.,11./48.,11./48.,1./36.},
		{23./48.,11./48., 3./48.,11./48.,1./36.},
		{23./48.,11./48.,11./48., 3./48.,1./36.},

		{11./48.,23./48., 3./48.,11./48.,1./36.},
		{11./48.,23./48.,11./48., 3./48.,1./36.},

		{11./48.,11./48.,23./48., 3./48.,1./36.},

		{ 3./48.,23./48.,11./48.,11./48.,1./36.},
		{ 3./48.,11./48.,23./48.,11./48.,1./36.},
		{ 3./48.,11./48.,11./48.,23./48.,1./36.},

		{11./48., 3./48.,23./48.,11./48.,1./36.},
		{11./48., 3./48.,11./48.,23./48.,1./36.},

		{11./48.,11./48., 3./48.,23./48.,1./36.},

		{ 7./48., 3./48.,19./48.,19./48.,1./36.},
		{ 7./48.,19./48., 3./48.,19./48.,1./36.},
		{ 7./48.,19./48.,19./48., 3./48.,1./36.},

		{19./48., 7./48., 3./48.,19./48.,1./36.},
		{19./48., 7./48.,19./48., 3./48.,1./36.},

		{19./48.,19./48., 7./48., 3./48.,1./36.},

		{ 3./48., 7./48.,19./48.,19./48.,1./36.},
		{ 3./48.,19./48., 7./48.,19./48.,1./36.},
		{ 3./48.,19./48.,19./48., 7./48.,1./36.},

		{19./48., 3./48., 7./48.,19./48.,1./36.},
		{19./48., 3./48.,19./48., 7./48.,1./36.},

		{19./48.,19./48., 3./48., 7./48.,1./36.}
	};
*/
	static struct{
		double coef[4];
		double weight;
	}volume_point[m_num_VP];

	static bool isfirst=true;
	if(isfirst){
		int cnt=0;
		for(int x=0;x<m_VP_div-3;x++)
		for(int y=0;y<m_VP_div-3-x;y++)
		for(int z=0;z<m_VP_div-3-x-y;z++){
			double dx=1./(double)(m_VP_div);
			volume_point[cnt].coef[0]=(double)(x+1)*dx;
			volume_point[cnt].coef[1]=(double)(y+1)*dx;
			volume_point[cnt].coef[2]=(double)(z+1)*dx;
			volume_point[cnt].coef[3]=(double)(m_VP_div-x-y-z-3)*dx;
			volume_point[cnt].weight=1./m_num_VP;
			cnt++;
		}
		isfirst=false;
	}

	if(!vfaces.size()){
		if(vmarker.size()){
			m_marker_Volume=m_Volume;
			m_weight=1.0;
		}else{
			m_marker_Volume=0;
			m_weight=0.0;
		}
		return;
	}

	double weight_sum=0;
	int cnt=0;
	for(int i=0;i<m_num_nodes;i++)
		m_shape_weight[i]=0;
	for(int vp=0;vp<m_num_VP;vp++){
		bool flag=true;
		double coordinate[3]={0.};
		for(int i=0;i<4;i++)
			for(int j=0;j<3;j++)
				coordinate[j]+=nodes[i]->dof[j].X*volume_point[vp].coef[i];

		for(unsigned int face=0;face<vfaces.size();face++)
			if(vfaces[face]->distance(coordinate)<0)
				flag=false;

		if(flag){
			weight_sum+=volume_point[vp].weight;
			for(int i=0;i<m_num_nodes;i++)
				m_shape_weight[i]+=volume_point[vp].coef[i];
			cnt++;
		}
	}
	for(int i=0;i<m_num_nodes;i++)
		m_shape_weight[i]*=(double)m_num_nodes/(double)cnt;

	//issue3
	if(!weight_sum)
		weight_sum=1./(double)m_num_VP;

#ifdef MAX_VOLUME
	if(0<weight_sum && weight_sum < 1.-1E-8)
		weight_sum=1.0;
#endif

	m_marker_Volume=m_Volume*weight_sum;
	m_weight=weight_sum;
}

void Tetra1_marker::getVolumeCoor(double dest[4],double x[3]){
	double p[3]={
		x[0]-nodes[3]->dof[0].X,
		x[1]-nodes[3]->dof[1].X,
		x[2]-nodes[3]->dof[2].X
	};
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			dest[i]+=m_dNdX[i][j]*p[j];
	dest[3]=1.-dest[0]-dest[1]-dest[2];
}

void Tetra1_marker::calcK_alpha(Sparse& K){
	for(int a=0;a<m_num_nodes;a++)
	for(int b=0;b<m_num_nodes;b++)
	for(int i=0;i<m_dim;i++){
		int row=nodes[a]->dof[i].mi;
		int col=nodes[b]->dof[i].mi;
		if(col<row || row<0)
			continue;
		double sum=0;
		for(int k=0;k<m_dim;k++)
			sum+=m_dNdx[a][k]*m_dNdx[b][k];
		K.add(row,col,
			//identity
			//m_diffuse*m_volume
			//heat transfer for each component
			sum*m_diffuse*m_Volume
			//heat transfer
			//m_dNdx[a][i]*m_dNdx[b][i]*m_diffuse*m_volume
		);
	}
}

void Tetra1_marker::calcK_c(Sparse& K){
	//fixme use iterator
	for(unsigned a=0;a<m_num_nodes;a++)
	for(unsigned b=0;b<m_num_nodes;b++)
	for(unsigned i=0;i<m_dim;i++)
	for(unsigned j=0;j<m_dim;j++){
		int row=nodes[a]->dof[i].mi;
		int col=nodes[b]->dof[j].mi;
		if(col<row || row<0)
			continue;
		if(!vmarker.size() || !m_marker_Volume){
			continue;
		}

		double numerator=0.0;
		double denominator=0.0;
		for(unsigned int m=0;m<vmarker.size();m++){
			denominator+=vmarker[m]->getWeight();
			for(int k=0;k<m_dim;k++)
			for(int l=0;l<m_dim;l++){
				numerator+=
					m_dNdx[a][k]*(
						vmarker[m]->c[i][k][j][l]
					)*m_dNdx[b][l]*vmarker[m]->getWeight();
			}
		}
		K.add(row,col,numerator/denominator*m_marker_volume);
	}
}

void Tetra1_marker::calcK_sigma(Sparse& K){
	for(unsigned a=0;a<m_num_nodes;a++)
	for(unsigned b=0;b<m_num_nodes;b++)
	for(unsigned i=0;i<m_dim;i++){
		int row=nodes[a]->dof[i].mi;
		int col=nodes[b]->dof[i].mi;
		if(col<row || row<0)
			continue;
		if(!vmarker.size() || !m_marker_Volume)
			continue;

		double denominator=0.0;
		double numerator=0.0;
		for(unsigned int m=0;m<vmarker.size();m++){
			denominator+=vmarker[m]->getWeight();
			for(int k=0;k<m_dim;k++)
			for(int l=0;l<m_dim;l++){
				numerator+=
					vmarker[m]->sigma(k,l)*
					m_dNdx[a][k]*
					m_dNdx[b][l]*
					vmarker[m]->getWeight();
			}
		}
		K.add(row,col,numerator/denominator*m_marker_volume);
	}
}

void Tetra1_marker::calcK_mass(Sparse& K){
	for(unsigned a=0;a<m_num_nodes;a++)
	for(unsigned b=0;b<m_num_nodes;b++)
	for(unsigned i=0;i<m_dim;i++){
		int row=nodes[a]->dof[i].mi;
		int col=nodes[b]->dof[i].mi;
		if(col<row || nodes[a]->dof[i].isDirichlet)
			continue;
		double coef=(a==b?2.0:1.0);
		K.add(row,col,coef/120.0*m_Volume*m_mass_coef);
	}
}

void Tetra1_marker::calcT_unique(MATRIX& T){
	for(int a=0;a<m_num_nodes;a++)
	for(int i=0;i<m_dim;i++){
		if(!vmarker.size() || !m_marker_Volume){
			continue;
		}

		double numerator=0.0;
		double denominator=0.0;
		for(unsigned int k=0;k<vmarker.size();k++){
			denominator+=vmarker[k]->getWeight();
			for(int j=0;j<m_dim;j++)
				numerator+=vmarker[k]->sigma(i,j)*m_dNdx[a][j]*vmarker[k]->getWeight();
		}
		T(nodes[a]->dof[i].fmi,0)+=numerator/denominator*m_marker_volume;
	}
#ifdef ALPHA
	for(int a=0;a<m_num_nodes;a++)
	for(int i=0;i<m_dim;i++){
		double sum=0;
		for(int b=0;b<m_num_nodes;b++)
		for(int j=0;j<m_dim;j++)
			sum+=m_dNdx[b][j]*nodes[b]->dof[i].u*m_dNdx[a][j];
		T(nodes[a]->dof[i].fmi,0)+=sum*m_Volume*m_diffuse;
	}
#endif
#ifdef MASS
	calcT_mass(T);
#endif
}

void Tetra1_marker::calcT_mass(MATRIX& T){
	for(unsigned a=0;a<m_num_nodes;a++)
	for(unsigned i=0;i<m_dim;i++){
		double sum=0;
		for(unsigned b=0;b<m_num_nodes;b++)
			sum+=nodes[b]->dof[i].u*(a==b?2.0:1.0);
		T(nodes[a]->dof[i].fmi,0)+=sum/120.0*m_Volume*m_mass_coef;
	}
}

