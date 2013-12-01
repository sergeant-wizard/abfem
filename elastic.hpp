/*
ver.03.03
file elastic.hpp
*/

#pragma once
#include "common.h"
#include "material.h"
#include "matrix_algebraic.h"

using ABFEM::MATH::kronecker;
using ABFEM::MATH::kronecker2;

class Venant_: public Elastic<3,MATRIX33>{
public:
	Venant_(const double E,const double nu):Elastic<3,MATRIX33>(E,nu){};
	void getC_e(double dest[][3][3][3],MATRIX33* F_in){
		for(unsigned i=0;i<3;i++)
		for(unsigned j=0;j<3;j++)
		for(unsigned k=0;k<3;k++)
		for(unsigned l=0;l<3;l++){
			dest[i][j][k][l]=
				lambda*kronecker2[i][j][k][l]
			+mu*(kronecker2[i][k][j][l]+kronecker2[i][l][j][k]);
		}
	};

	void getSigma(MATRIX33* dest,MATRIX33* F_in){
		Matrix33 Fd(F_in);
		double J=Fd.I1()+Fd.I2()+Fd.det()+1.;
		Matrix33 bd=Fd.transpose()+Fd+Fd*Fd.transpose();
		Matrix33 bd2=bd*bd;

		for(unsigned i=0;i<3;i++)
		for(unsigned j=0;j<3;j++)
			(*dest)(i,j)=
			0.5*lambda/J*bd.I1()
			*(bd(i,j)+kronecker[i][j])
			+mu/J*(bd2(i,j)+bd(i,j));
	};
};

class Bonet_comp_: public Elastic<3,MATRIX33>{
public:
	Bonet_comp_(const double E,const double nu):Elastic<3,MATRIX33>(E,nu){};
	void getC_e(double dest[][3][3][3],MATRIX33* F_in){
		Matrix33 Fd(F_in);
		double J_=Fd.I1()+Fd.I2()+Fd.det();
		for(unsigned i=0;i<3;i++)
		for(unsigned j=0;j<3;j++)
		for(unsigned k=0;k<3;k++)
		for(unsigned l=0;l<3;l++){
			dest[i][j][k][l]=(
				lambda*kronecker2[i][j][k][l]+
				(mu-lambda*log1p(J_))*(kronecker2[i][k][j][l]+kronecker2[i][l][j][k])
			)/(1.+J_);
		}
	};

	void getSigma(MATRIX33* dest,MATRIX33* F_in){
		Matrix33 Fd(F_in);
		double J_=Fd.I1()+Fd.I2()+Fd.det();
		Matrix33 bd=Fd.transpose()+Fd+Fd*Fd.transpose();

		for(unsigned i=0;i<3;i++)
		for(unsigned j=0;j<3;j++)
			(*dest)(i,j)=(
				mu*bd(i,j)+lambda*log1p(J_)*kronecker[i][j]
			)/(1.+J_);
	};
};
