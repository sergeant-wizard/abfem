/*
ver.03.03
file incompressible.hpp
*/

#pragma once
#include <math.h>
#include "material.h"
#include "common.h"
#include "matrix_algebraic.h"

using ABFEM::MATH::kronecker;
using ABFEM::MATH::kronecker2;

class Neo_hookean:public Incompressible<3,MATRIX33>{
public:
	Neo_hookean(const double E):
		Elastic<3,MATRIX33>(E,0.5),
		Incompressible<3,MATRIX33>(E)
	{};
	void getC_e(double dest[][3][3][3],MATRIX33* F_in){
		Matrix33 F(F_in);
		for(unsigned i=0;i<3;i++)
			F(i,i)+=1.0;
		double J=fabs(F.I3());
		Matrix33 b=F*F.transpose();
		double I_b=b.I1();
		for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
		for(int k=0;k<3;k++)
		for(int l=0;l<3;l++){
			dest[i][j][k][l]=
				2.0*mu*pow(J,-5.0/3.0)*(//sigma-base
				//2.0*mu*pow(J,-2.0/3.0)*(//tau-base
					I_b*(kronecker2[i][k][j][l]+kronecker2[i][l][j][k])/6.0
					-b(i,j)*kronecker[k][l]/3.0
					-b(k,l)*kronecker[i][j]/3.0
					+I_b*kronecker2[i][j][k][l]/9.0
				);
		}
	};

	void getSigma(MATRIX33* dest,MATRIX33* F_in){
		Matrix33 F_(F_in);
		double J=F_.I1()+F_.I2()+F_.I3()+1.;
		Matrix33 b=F_*F_.transpose()+F_+F_.transpose();
		double I_b=b.I1();
		for(int i=0;i<3;i++)
			for(int j=0;j<3;j++)
				(*dest)(i,j)=mu*pow(J,-5.0/3.0)* //sigma-base
				//(*dest)(i,j)=mu*pow(J,-2.0/3.0)* //tau-base
					(b(i,j)-I_b*kronecker[i][j]/3.0);
	};

	//double potential(const MATRIX33* F_)const{
	double potential(const MATRIX33* F_in){
		using ABFEM::MATH::SQR;
		Matrix33 F_(F_in);
		double J=F_.I1()+F_.I2()+F_.I3()+1.0;
		double detC3=pow(J,2./3.);
		double ret=0;
		for(unsigned i=0;i<3;i++)
		for(unsigned j=0;j<3;j++)
			ret+=SQR(F_(i,j));
		for(unsigned i=0;i<3;i++)
			ret+=2.0*F_(i,i);
		ret+=3.0*(1.0-detC3);
		ret*=0.5*mu/detC3;
		if(isnan(ret))
			ret=-1;
		return ret;

	};
};

class Mooney_Rivlin:public Incompressible<3,MATRIX33>{
public:
	Mooney_Rivlin(const double c1,const double c2):
		Elastic<3,MATRIX33>((c1+c2)*6.,0.5),
		Incompressible<3,MATRIX33>((c1+c2)*6.),
		c1(c1),
		c2(c2)
	{};
	void getC_e(double dest[][3][3][3],MATRIX33* F_in){
		Matrix33 F(F_in);
		for(unsigned i=0;i<3;i++)
			F(i,i)+=1.0;
		double J=fabs(F.I3());

		Matrix33 C=F.transpose()*F;
		double I_C  =C.I1();
		double II_C =C.I2();
		double III_C=C.I3();

		Matrix33 C_inv=C.inv();

		double W1=c1*pow(III_C,-1./3.);
		double W2=c2*pow(III_C,-2./3.);

		const static double four_ninth=4./9.;
		double c_material[3][3][3][3];

		for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
		for(int k=0;k<3;k++)
		for(int l=0;l<3;l++){
			c_material[i][j][k][l]=four_ninth*(
				-9.*W2*.5*(kronecker2[i][k][j][l]+kronecker2[i][l][j][k])	//I4
				+9.*W2*kronecker2[i][j][k][l]								//I * I
				+(W1*I_C+4.*W2*II_C)*C_inv(i,j)*C_inv(k,l)					//C-1*C-1
				+3.*(W1*I_C+2.*W2*II_C)*.5*(								//C/C
					C_inv(i,k)*C_inv(j,l)+
					C_inv(i,l)*C_inv(j,k)
				)
				-3.*(W1+2.*W2*I_C)*(										//C-1*I
					C_inv(i,j)*kronecker[k][l]+
					C_inv(k,l)*kronecker[i][j]
				)
				+6.*W2*(
					C(i,j)*C_inv(k,l)*
					C(k,l)*C_inv(i,j)
				)
			);
		}

		for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
		for(int k=0;k<3;k++)
		for(int l=0;l<3;l++){
			double sum=0;
			for(int I=0;I<3;I++)
			for(int J=0;J<3;J++)
			for(int K=0;K<3;K++)
			for(int L=0;L<3;L++)
				sum+=F(i,I)*F(j,J)*F(k,K)*F(l,L)*c_material[I][J][K][L]/J;
			dest[i][j][k][l]=sum;
		}
	};

	void getSigma(MATRIX33* dest,MATRIX33* F_in){
		Matrix33 F(F_in);
		for(unsigned i=0;i<3;i++)
			F(i,i)+=1.0;
		double J=F.I3();

		Matrix33 F_T=F.transpose();
		Matrix33 C=F_T*F;
		double I_C  =C.I1();
		double II_C =C.I2();
		double III_C=C.I3();

		Matrix33 C_inv=C.inv();

		Matrix33 I;
		I.identity();

		Matrix33 c1_mat,c2_mat;
		{
			c1_mat=-1./3.*I_C*C_inv;
			for(unsigned i=0;i<3;i++)
				c1_mat(i,i)+=1.;
			c1_mat*=c1*pow(III_C,-1./3.);
		}{
			c2_mat=I_C*I-C-2./3.*II_C*C_inv;
			c2_mat*=c2*pow(III_C,-2./3.);
		}

		Matrix33 S=c1_mat+c2_mat;
		Matrix33 tau=F*S*F_T;

		for(int i=0;i<3;i++)
			for(int j=0;j<3;j++)
				(*dest)(i,j)=tau(i,j)/J*2.;
	};

	double potential(const MATRIX33* F_in){
		Matrix33 F(F_in);
		for(unsigned i=0;i<3;i++)
			F(i,i)+=1.0;

		Matrix33 C=F.transpose()*F;
		double I_C  =F.I1();
		double II_C =F.I2();
		double III_C=F.I3();
		return
			c1*F.I1()*pow(III_C,-1./3.)+
			c2*F.I2()*pow(III_C,-2./3.);
	};
private:
	const double c1;
	const double c2;
};
