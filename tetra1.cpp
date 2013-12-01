/*
ver.03.03
file tetra1.cpp
*/

#include "tetra1.h"
using std::vector;

Tetra1::Tetra1(int* node_index,Node* raw_nodes,Default_material* material):
	material(material),
	m_Volume(0)
{
	m_F.identity();
	m_F_prev.identity();
	m_F_rel.identity();
	m_F_.zero_out();
	m_F_rel_.zero_out();

	for(int i=0;i<m_num_nodes;i++)
		nodes[i]=&raw_nodes[node_index[i]];
}

Tetra1::~Tetra1(){
}

void Tetra1::initialize(void){
	update_loadstep();
}

void Tetra1::initialize_sparse(Sparse& K){
	for(int a=0;a<m_num_nodes;a++)
		for(int b=0;b<m_num_nodes;b++)
			for(int i=0;i<m_dim;i++){
				int row=nodes[a]->dof[i].mi;
				if(nodes[a]->dof[i].isDirichlet)
					continue;
				for(int j=0;j<m_dim;j++){
					int col=nodes[b]->dof[j].mi;
					if(col<row)
						continue;
					K.add(row,col);
				}
			}
}

void Tetra1::update_loadstep(void){
	double& x1=nodes[0]->dof[0].X;
	double& y1=nodes[0]->dof[1].X;
	double& z1=nodes[0]->dof[2].X;
	double& x2=nodes[1]->dof[0].X;
	double& y2=nodes[1]->dof[1].X;
	double& z2=nodes[1]->dof[2].X;
	double& x3=nodes[2]->dof[0].X;
	double& y3=nodes[2]->dof[1].X;
	double& z3=nodes[2]->dof[2].X;
	double& x4=nodes[3]->dof[0].X;
	double& y4=nodes[3]->dof[1].X;
	double& z4=nodes[3]->dof[2].X;

	double Det=
		x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + 
		x2*y3*z1 + x3*y1*z2 - x3*y2*z1 - 
		x1*y2*z4 + x1*y4*z2 + x2*y1*z4 - 
		x2*y4*z1 - x4*y1*z2 + x4*y2*z1 + 
		x1*y3*z4 - x1*y4*z3 - x3*y1*z4 + 
		x3*y4*z1 + x4*y1*z3 - x4*y3*z1 - 
		x2*y3*z4 + x2*y4*z3 + x3*y2*z4 - 
		x3*y4*z2 - x4*y2*z3 + x4*y3*z2;

	//issue10
#ifdef TOTAL_LAGRANGE
	if(!m_Volume )
		m_Volume=fabs(Det)/6.;
#else
	m_Volume=fabs(Det)/6.;
#endif

	m_dNdX[0][0]=(y2*z3 - y3*z2 - y2*z4 + y4*z2 + y3*z4 - y4*z3)/Det;
	m_dNdX[0][1]=(x3*z2 - x2*z3 + x2*z4 - x4*z2 - x3*z4 + x4*z3)/Det;
	m_dNdX[0][2]=(x2*y3 - x3*y2 - x2*y4 + x4*y2 + x3*y4 - x4*y3)/Det;
	m_dNdX[1][0]=(y3*z1 - y1*z3 + y1*z4 - y4*z1 - y3*z4 + y4*z3)/Det;
	m_dNdX[1][1]=(x1*z3 - x3*z1 - x1*z4 + x4*z1 + x3*z4 - x4*z3)/Det;
	m_dNdX[1][2]=(x3*y1 - x1*y3 + x1*y4 - x4*y1 - x3*y4 + x4*y3)/Det;
	m_dNdX[2][0]=(y1*z2 - y2*z1 - y1*z4 + y4*z1 + y2*z4 - y4*z2)/Det;
	m_dNdX[2][1]=(x2*z1 - x1*z2 + x1*z4 - x4*z1 - x2*z4 + x4*z2)/Det;
	m_dNdX[2][2]=(x1*y2 - x2*y1 - x1*y4 + x4*y1 + x2*y4 - x4*y2)/Det;
	m_dNdX[3][0]=(y2*z1 - y1*z2 + y1*z3 - y3*z1 - y2*z3 + y3*z2)/Det;
	m_dNdX[3][1]=(x1*z2 - x2*z1 - x1*z3 + x3*z1 + x2*z3 - x3*z2)/Det;
	m_dNdX[3][2]=(x2*y1 - x1*y2 + x1*y3 - x3*y1 - x2*y3 + x3*y2)/Det;

	m_F_prev.copy(m_F);
	m_F_prev_.copy(m_F_);

	m_sigma_prev.copy(m_sigma);
}

void Tetra1::update_F(void){
	//update deformation gradient
	for(int i=0;i<m_dim;i++){
		for(int j=0;j<m_dim;j++){
			double sum=0;
			for(int k=0;k<m_num_nodes;k++)
				sum+=nodes[k]->dof[i].x*m_dNdX[k][j];
			m_F_rel(i,j)=sum;
		}
	}
	m_F=m_F_rel*m_F_prev;
	m_J=m_F.I3();

	//update deformation gradient_
	for(int i=0;i<m_dim;i++){
		for(int j=0;j<m_dim;j++){
			double sum=0;
			for(int k=0;k<m_num_nodes;k++)
				sum+=nodes[k]->dof[i].u*m_dNdX[k][j];
			m_F_rel_(i,j)=sum;
		}
	}
	m_F_=m_F_rel_*m_F_prev_+m_F_rel_+m_F_prev_;
	m_J_=m_F_.I3()+m_F_.I2()+m_F_.I1();
	m_J_rel_=m_F_rel_.I3()+m_F_rel_.I2()+m_F_rel_.I1();

	m_volume=m_Volume*(1.+m_J_rel_);
};

void Tetra1::update_Constitutive(void){
	material->getC_e(m_c,&m_F_);
	material->getSigma(&m_sigma,&m_F_);
}

void Tetra1::update_dNdx(void){
	Matrix33 Finv=m_F_rel.inv();
	double* pdNdx=m_dNdx[0];
	for(int node=0;node<m_num_nodes;node++){
		for(int component=0;component<m_dim;component++){
			*pdNdx=0.0;
			for(int i=0;i<m_dim;i++)
				*pdNdx+=Finv(i,component)*m_dNdX[node][i];
			pdNdx++;
		}
	}
}

void Tetra1::update_iteration(void){
	update_F();
	update_Constitutive();
	update_dNdx();
}

void Tetra1::sync_Kmatrix(Sparse& K,bool symmetric){
	vmatrix_iterator.clear();
	for(int a=0;a<m_num_nodes;a++)
	for(int b=0;b<m_num_nodes;b++)
	for(int i=0;i<m_dim;i++)
	for(int j=0;j<m_dim;j++){
		int row=nodes[a]->dof[i].mi;
		int col=nodes[b]->dof[j].mi;
		if( (symmetric && col<row) || nodes[a]->dof[i].isDirichlet )
			continue;
		vmatrix_iterator.push_back(
			Matrix_iterator(K(row,col),a,b,i,j)
		);
	}
}

void Tetra1::calcK_c(Sparse& K){
#ifdef ITERATOR
	for(vector<Matrix_iterator>::iterator it=vmatrix_iterator.begin();it!=vmatrix_iterator.end();it++){
		double sum=0.0;
		for(int k=0;k<m_dim;k++)
			for(int l=0;l<m_dim;l++)
				sum+=m_dNdx[it->a][k]*m_dNdx[it->b][l]*m_c[it->i][k][it->j][l];
		it->ref+=sum*m_volume;
	}
#else
	for(int a=0;a<m_num_nodes;a++)
	for(int b=0;b<m_num_nodes;b++)
	for(int i=0;i<m_dim;i++)
	for(int j=0;j<m_dim;j++){
		MATRIXINDEX row=nodes[a]->dof[i].mi;
		MATRIXINDEX col=nodes[b]->dof[j].mi;
		if(col<row || nodes[a]->dof[i].isDirichlet)
			continue;
		double sum=0.0;
		for(int k=0;k<m_dim;k++)
			for(int l=0;l<m_dim;l++)
				sum+=m_dNdx[a][k]*m_dNdx[b][l]*m_c[i][k][j][l];
		K.add(row,col,sum*m_volume);
	}
#endif
	//wth mass
#ifdef MASS
	using ABFEM::CONSTANTS::g_mass_coef;
	static const double m_mass_coef=g_mass_coef;
	for(unsigned a=0;a<m_num_nodes;a++)
	for(unsigned b=0;b<m_num_nodes;b++)
	for(unsigned i=0;i<m_dim;i++){
		int row=nodes[a]->dof[i].mi;
		int col=nodes[b]->dof[i].mi;
		if(col<row || nodes[a]->dof[i].isDirichlet)
			continue;
		double coef=(a==b?2.0:1.0);
		K.add(row,col,coef/120.0*m_volume*m_mass_coef);
	}
#endif
}

void Tetra1::calcK_sigma(Sparse& K){
	for(vector<Matrix_iterator>::iterator it=vmatrix_iterator.begin();it!=vmatrix_iterator.end();it++){
		if(it->i != it->j)
			continue;
		double sum=0;
		for(int k=0;k<m_dim;k++)
			for(int l=0;l<m_dim;l++)
				sum+=
				m_sigma(k,l)*
				m_dNdx[it->a][k]*
				m_dNdx[it->b][l];
		it->ref+=sum*m_volume;
	}
}

void Tetra1::calcT(MATRIX& T){
	for(int a=0;a<m_num_nodes;a++)
	for(int i=0;i<m_dim;i++){
		double sum=0;
		for(int j=0;j<m_dim;j++){
			sum+=m_sigma(i,j)*m_dNdx[a][j];
		}
		T(nodes[a]->dof[i].fmi,0)+=sum*m_volume;
	}

	//wth mass
#ifdef MASS
	using ABFEM::CONSTANTS::g_mass_coef;
	static const double m_mass_coef=g_mass_coef;
	for(unsigned a=0;a<m_num_nodes;a++)
	for(unsigned i=0;i<m_dim;i++){
		double sum=0;
		for(unsigned b=0;b<m_num_nodes;b++)
			sum+=nodes[b]->dof[i].u*(a==b?2.0:1.0);
		T(nodes[a]->dof[i].fmi,0)+=sum/120.0*m_volume*m_mass_coef;
	}
#endif
}

//issue12
void Tetra1::calcBodyForce(MATRIX& T){
	extern double g_density;
	extern double g_gravity;
	extern double g_bodyforce_normal[3];
	extern int g_loadstep_cnt;
	extern int g_loadsteps;
	for(int a=0;a<m_num_nodes;a++)
	for(int i=0;i<m_dim;i++){
		T(nodes[a]->dof[i].fmi,0)+=g_density*g_gravity*m_Volume*0.25*g_bodyforce_normal[i]*((double)g_loadstep_cnt/(double)g_loadsteps);
	}
}

void Tetra1::debug_(void){
/*
	for(int i=0;i<4;i++){
		for(int j=0;j<3;j++){
			cout << nodes[i]->ci[j].x << "\t";
		}
		cout<<endl;
	}

	for(int i=0;i<4;i++){
		for(int j=0;j<3;j++){
			cout << nodes[i]->ci[j].X << "\t";
		}
		cout<<endl;
	}

	cout << "elem volume" << endl;
	cout << (m_volume) << endl;
	cout << "volume" << endl;
	cout << m_Volume << endl;

*/
/*
	for(int i=0;i<vfaces.size();i++){
		cout<<"face"<<i<<endl;
		double* vertice;
		for(int j=0;j<3;j++){
			vertice=vfaces[i]->getVertice(j);
			for(int k=0;k<3;k++)
				cout << vertice[k] <<"\t";
			cout<<endl;
		}
	}
		
	cout<<vmarker.size() <<"\t"<<vfaces.size();
*/
	//(material.*(material.getSigma))(sigma,m_F);
	//cout<<m_F<<endl;
	//cout<<sigma(2,2)<<endl;
	//cout<<sigma<<endl;

	//cout << m_F << endl;
/*
	for(int i=0;i<4;i++){
		for(int j=0;j<3;j++){
			cout << this->nodes[i]->ci[j].X << "\t";
		}
		cout<<endl;
	}
*/
/*
	MATRIX33 K1,K2;
	double c_[3][3][3][3];
	(material.*(material.getC_e))(c_,m_F);
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			double sum=0;
			for(int k=0;k<3;k++){
				for(int l=0;l<3;l++){
					sum+=c_[i][k][j][l]*
					(m_dNdx[0][k]+kronecker[i][k]*(
						parent->get_dndx(super_index[0],k)-m_dNdx[0][k]
					)/3.)*
					(m_dNdx[1][l]+kronecker[j][l]*(
						parent->get_dndx(super_index[1],l)-m_dNdx[1][l]
					)/3.);
				}
			}
			K1(i,j)=sum;
		}
	}

	MATRIX B0(6,3);
	B0.Input(18,
		m_dNdx[0][0]+(parent->get_dndx(super_index[0],0)-m_dNdx[0][0])/3., 0., 0.,
		0., m_dNdx[0][1]+(parent->get_dndx(super_index[0],1)-m_dNdx[0][1])/3., 0.,
		0., 0., m_dNdx[0][2]+(parent->get_dndx(super_index[0],2)-m_dNdx[0][2])/3.,
		m_dNdx[0][1],m_dNdx[0][0],0.,
		m_dNdx[0][2],0.,m_dNdx[0][0],
		0.,m_dNdx[0][2],m_dNdx[0][1]
	);
	MATRIX B1(6,3);
	B1.Input(18,
		m_dNdx[1][0]+(parent->get_dndx(super_index[1],0)-m_dNdx[1][0])/3., 0., 0.,
		0., m_dNdx[1][1]+(parent->get_dndx(super_index[1],1)-m_dNdx[1][1])/3., 0.,
		0., 0., m_dNdx[1][2]+(parent->get_dndx(super_index[1],2)-m_dNdx[1][2])/3.,
		m_dNdx[1][1],m_dNdx[1][0],0.,
		m_dNdx[1][2],0.,m_dNdx[1][0],
		0.,m_dNdx[1][2],m_dNdx[1][1]
	);
*/
	/*
	MATRIX D(6,6);
	D.Input(36,
		m_c[0][0][0][0],m_c[0][0][1][1],m_c[0][0][2][2],m_c[0][0][0][1],m_c[0][0][0][2],m_c[0][0][1][2],
		m_c[1][1][0][0],m_c[1][1][1][1],m_c[1][1][2][2],m_c[1][1][0][1],m_c[1][1][0][2],m_c[1][1][1][2],
		m_c[2][2][0][0],m_c[2][2][1][1],m_c[2][2][2][2],m_c[2][2][0][1],m_c[2][2][0][2],m_c[2][2][1][2],
		m_c[0][1][0][0],m_c[0][1][1][1],m_c[0][1][2][2],m_c[0][1][0][1],m_c[0][1][0][2],m_c[0][1][1][2],
		m_c[0][2][0][0],m_c[0][2][1][1],m_c[0][2][2][2],m_c[0][2][0][1],m_c[0][2][0][2],m_c[0][2][1][2],
		m_c[1][2][0][0],m_c[1][2][1][1],m_c[1][2][2][2],m_c[1][2][0][1],m_c[1][2][0][2],m_c[1][2][1][2]
	);

	static bool isFirst=true;
	if(isFirst){
		isFirst=false;
		cout<<D<<endl;
	}
	*/
		

/*
	MATRIX tmp(3,6);
	tmp=B0.Transpose()*D;
*/
}

