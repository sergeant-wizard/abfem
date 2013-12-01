/*
ver.03.03
file matrix_algebraic.h
*/

#pragma once
#include <algorithm>
#include <iostream>
#include "matrix.h"

class Matrix_body{
	typedef std::ostream ostream;
	friend class Matrix;
public:
	virtual ~Matrix_body();
private:
	Matrix_body(unsigned row,unsigned col);
	double& operator()(unsigned i,unsigned j);
	double operator()(unsigned i,unsigned j)const;
	void print(ostream& out)const;
	short referenceCount;
	const unsigned row;
	const unsigned col;
	double* val;
};

class Matrix:virtual public MATRIX{
	typedef std::ostream ostream;
public:
	Matrix(void);
	Matrix(const Matrix& src);
	Matrix(unsigned row,unsigned col);
	virtual ~Matrix();
	unsigned GetRow(void)const;
	unsigned GetCol(void)const;
	Matrix& operator=(const Matrix& right);
	double& operator()(unsigned i,unsigned j);
	double operator()(unsigned i,unsigned j)const;
	double& operator[](unsigned i);
	Matrix operator-(void)const;
	Matrix operator+(Matrix right)const;
	Matrix operator-(Matrix right)const;
	Matrix operator*(Matrix right)const;
	Matrix operator*(double coef)const;
	Matrix& operator*=(double coef);
	friend Matrix operator*(double coef,Matrix right);
	Matrix transpose(void)const;

	void copy(const Matrix& src);
	void zero_out(void);
	void print(ostream& out)const;
protected:
	Matrix_body* rep;
};

class Matrix33:public Matrix,public MATRIX33{
public:
	Matrix33(void):Matrix(3,3){};
	Matrix33(const Matrix& src):Matrix(src){};
	Matrix33(const MATRIX33* src):Matrix(3,3){
		for(unsigned i=0;i<3;i++)
		for(unsigned j=0;j<3;j++)
			(*this)(i,j)=(*src)(i,j);
	};
	Matrix& operator=(const Matrix& right){return Matrix::operator=(right);};

	double I1(void)const;
	double I2(void)const;
	double I3(void)const;
	void identity(void);
	void dev(void);
	Matrix33 inv(void)const;
};

