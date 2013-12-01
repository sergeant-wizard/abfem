/*
ver.03.03
file matrix.h
*/

#pragma once
#include <string>

class MATRIX{
public:
	MATRIX(void){};
	MATRIX(MATRIX*){};
	virtual ~MATRIX(){};

	//virtual void Renew(unsigned row,unsigned col)=0;
	//virtual void ZeroOut(void)=0;
	//virtual void Identity(void)=0;
	//virtual double* GetData(void)=0;
	//virtual double* GetData(void)const=0;
	//virtual void Print(std::string filename)=0;

	virtual unsigned GetRow(void)const=0;
	virtual unsigned GetCol(void)const=0;

	virtual double& operator()(unsigned i,unsigned j)=0;
	virtual double operator()(unsigned i,unsigned j)const=0;
	virtual double& operator[](unsigned i)=0;
	MATRIX& operator=(MATRIX& right){return *this;};
};

class MATRIX33:virtual public MATRIX{
public:
	MATRIX33(void){};
	MATRIX33(MATRIX* src){};
	virtual ~MATRIX33(){};
	MATRIX33& operator=(MATRIX33& right){
		for(unsigned i=0;i<9;i++)
			(*this)[i]=right[i];
		return *this;
	};
	virtual double det(void)const{
		return
			(*this)(0,0)*(*this)(1,1)*(*this)(2,2) +
			(*this)(0,1)*(*this)(1,2)*(*this)(2,0) +
			(*this)(0,2)*(*this)(1,0)*(*this)(2,1) -
			(*this)(0,2)*(*this)(1,1)*(*this)(2,0) -
			(*this)(0,0)*(*this)(1,2)*(*this)(2,1) -
			(*this)(0,1)*(*this)(1,0)*(*this)(2,2);
	};
	double I1(void)const{
		return (*this)(0,0)+(*this)(1,1)+(*this)(2,2);
	};
	double I2(void)const{
		double sum=0;
		for(unsigned i=0;i<3;i++)
		for(unsigned j=0;j<3;j++)
			sum+=(*this)(i,j)*(*this)(j,i);
		double I1=this->I1();
		return (I1*I1-sum)*.5;
	};
	void Inv(MATRIX33& dest){
		double mdet=det();
		dest[0]=((*this)(1,1)*(*this)(2,2) - (*this)(1,2)*(*this)(2,1))/mdet;
		dest[1]=((*this)(0,2)*(*this)(2,1) - (*this)(0,1)*(*this)(2,2))/mdet;
		dest[2]=((*this)(0,1)*(*this)(1,2) - (*this)(0,2)*(*this)(1,1))/mdet;
		dest[3]=((*this)(1,2)*(*this)(2,0) - (*this)(1,0)*(*this)(2,2))/mdet;
		dest[4]=((*this)(0,0)*(*this)(2,2) - (*this)(0,2)*(*this)(2,0))/mdet;
		dest[5]=((*this)(0,2)*(*this)(1,0) - (*this)(0,0)*(*this)(1,2))/mdet;
		dest[6]=((*this)(1,0)*(*this)(2,1) - (*this)(1,1)*(*this)(2,0))/mdet;
		dest[7]=((*this)(0,1)*(*this)(2,0) - (*this)(0,0)*(*this)(2,1))/mdet;
		dest[8]=((*this)(0,0)*(*this)(1,1) - (*this)(0,1)*(*this)(1,0))/mdet;
	};
private:
};

