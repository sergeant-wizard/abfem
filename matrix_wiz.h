/*
ver.03.03
file matrix_wiz.h
*/

#pragma once
#include <cstdarg>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include "matrix.h"

class MATRIX_WIZ:virtual public MATRIX{
public:
	MATRIX_WIZ();
	MATRIX_WIZ(unsigned row,unsigned col);
	MATRIX_WIZ(MATRIX_WIZ& src);
	MATRIX_WIZ(MATRIX_WIZ* src);
	virtual ~MATRIX_WIZ();
	void Renew(unsigned row,unsigned col);

	//utilities
	void ZeroOut(void);
	bool IsSqr(void);
	void Input(unsigned num,...);
	void Divide(MATRIX_WIZ& smaller,unsigned top,unsigned bottom,unsigned left, unsigned right);
	void Identity(void);
	void Print(std::string filename);
	unsigned GetRow(void)const{return row;}
	unsigned GetCol(void)const{return col;}

	//calculation
	void Eigen(MATRIX_WIZ& eigenvectors,MATRIX_WIZ& eigenvalues);
	double Trace(void);
	double Norm2(void);
	void Inv(MATRIX_WIZ& dest);
	MATRIX_WIZ& Transpose(void);

	//operators
	MATRIX& operator=(MATRIX_WIZ& right);
	MATRIX& operator+=(MATRIX& right);
	MATRIX& operator-=(MATRIX& right);
	MATRIX& operator*=(double coef);
	MATRIX& operator*=(MATRIX& right);
	friend std::ostream& operator<<(std::ostream& os,MATRIX_WIZ& matrix);
	friend MATRIX_WIZ& operator+(MATRIX_WIZ& left,MATRIX_WIZ& right);
	friend MATRIX_WIZ& operator-(MATRIX_WIZ& left,MATRIX_WIZ& right);
	friend MATRIX_WIZ& operator*(MATRIX_WIZ& mat1,MATRIX_WIZ& mat2);
	double& operator()(unsigned i,unsigned j);
	double operator()(unsigned i,unsigned j)const;
	double& operator[](unsigned i){return data[i];};
	double Atm(unsigned rown,unsigned coln);
	double* GetData(void){return &data[0];};
	double* GetData(void)const{return &data[0];};
#ifndef _DEBUG
protected:
#endif
	unsigned index(unsigned _row,unsigned _col);
	void gauss(unsigned _row,double coef);
	void gauss(unsigned row1,unsigned row2,double coef);
	void gauss(unsigned row1,unsigned row2);
	unsigned row;
	unsigned col;
	unsigned size;
	double* data;
	double* data_cnt;
	bool input_flag;
	const static unsigned width=8;

	static const double m_eps;
};

class MATRIX33_WIZ:public MATRIX_WIZ,public MATRIX33{
public:
	MATRIX33_WIZ();
	MATRIX33_WIZ(MATRIX33&);
	MATRIX33_WIZ(MATRIX33*);
	~MATRIX33_WIZ(){
		if(data)
			delete[] data;
		data=NULL;
	};
	void Inv(MATRIX33& dest);

	MATRIX& operator=(MATRIX33_WIZ& right);
	MATRIX& operator=(MATRIX_WIZ& right);
	MATRIX33_WIZ& Transpose(void);
	void dev(void);

	/*obsolete
	double det(void)const;
	double I1(void)const;
	double I2(void)const;
	*/
};

