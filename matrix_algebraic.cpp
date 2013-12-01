/*
ver.03.03
file matrix_algebraic.cpp
*/

#include "matrix_algebraic.h"

#ifdef _DEBUG
using namespace std;
#endif

Matrix_body::Matrix_body(unsigned row,unsigned col):
	referenceCount(1),
	row(row),col(col),
	val(0)
{
#ifdef _DEBUG
	cout << "body construtor called" << endl;
#endif
	val=new double[row*col];
}

Matrix_body::~Matrix_body(){
#ifdef _DEBUG
	cout << "body destructor called" << endl;
#endif
	if(val)
		delete[] val;
}

double& Matrix_body::operator()(unsigned i,unsigned j){
	return val[i*col+j];
}

double Matrix_body::operator()(unsigned i,unsigned j)const{
	return val[i*col+j];
}

void Matrix_body::print(ostream& out)const{
	for(unsigned i=0;i<row;i++){
		for(unsigned j=0;j<col;j++)
			out << (*this)(i,j) << " ";
		out<<std::endl;
	}
}

Matrix::Matrix(void):
	rep(0)
{
}

Matrix::Matrix(const Matrix& src):
	rep(src.rep)
{
	if(src.rep)
		src.rep->referenceCount++;
}

Matrix::Matrix(unsigned row,unsigned col):
	rep(0)
{
	rep=new Matrix_body(row,col);
}

Matrix::~Matrix(){
	if(rep && --rep->referenceCount==0)
		delete rep;
}

unsigned Matrix::GetRow(void)const{
	if(rep)
		return rep->row;
	else
		return 0;
}

unsigned Matrix::GetCol(void)const{
	if(rep)
		return rep->col;
	else
		return 0;
}

Matrix& Matrix::operator=(const Matrix& right){
	Matrix temp(right);
	std::swap(this->rep,temp.rep);
	return *this;
}

double& Matrix::operator()(unsigned i,unsigned j){
	static double dummy;
	if(rep)
		return (*rep)(i,j);
	else
		return dummy;
}

double Matrix::operator()(unsigned i,unsigned j)const{
	if(rep)
		return (*rep)(i,j);
	else
		return 0;
}

double& Matrix::operator[](unsigned i){
	static double dummy;
	if(rep)
		return (*rep)(i/rep->col,i%rep->col);
	else
		return dummy;
}

void Matrix::copy(const Matrix& src){
	if(!src.rep)
		return;
	Matrix temp(src.GetRow(),src.GetCol());
	for(unsigned i=0;i<src.GetRow();i++)
		for(unsigned j=0;j<src.GetRow();j++)
			temp(i,j)=src(i,j);
	std::swap(*this,temp);
}

void Matrix::zero_out(void){
	for(unsigned i=0;i<rep->row;i++)
	for(unsigned j=0;j<rep->col;j++)
		(*this)(i,j)=0;
}

void Matrix::print(ostream& out)const{
	rep->print(out);
}

Matrix Matrix::operator-(void)const{
	Matrix temp(rep->row,rep->col);
	for(unsigned i=0;i<rep->row;i++)
		for(unsigned j=0;j<rep->col;j++)
			temp(i,j)=-(*this)(i,j);
	return temp;
}

Matrix Matrix::operator+(Matrix right)const{
	Matrix temp(rep->row,rep->col);
	for(unsigned i=0;i<rep->row;i++)
		for(unsigned j=0;j<rep->col;j++)
			temp(i,j)=(*this)(i,j)+right(i,j);
	return temp;
}

Matrix Matrix::operator-(Matrix right)const{
	Matrix temp(rep->row,rep->col);
	for(unsigned i=0;i<rep->row;i++)
		for(unsigned j=0;j<rep->col;j++)
			temp(i,j)=(*this)(i,j)-right(i,j);
	return temp;
}

Matrix Matrix::operator*(Matrix right)const{
	Matrix temp(rep->row,right.rep->col);
	for(unsigned i=0;i<rep->row;i++)
		for(unsigned j=0;j<right.rep->col;j++){
			double sum=0;
			for(unsigned k=0;k<rep->col;k++)
				sum+=(*this)(i,k)*right(k,j);
			temp(i,j)=sum;
		}
	return temp;
}

Matrix Matrix::operator*(double coef)const{
	Matrix temp(rep->row,rep->col);
	for(unsigned i=0;i<rep->row;i++)
		for(unsigned j=0;j<rep->col;j++)
			temp(i,j)=(*this)(i,j)*coef;
	return temp;
}

Matrix& Matrix::operator*=(double coef){
	for(unsigned i=0;i<rep->row;i++)
		for(unsigned j=0;j<rep->col;j++)
			(*this)(i,j)*=coef;
	return *this;
}

Matrix operator*(double coef,Matrix right){
	Matrix temp(right.GetRow(),right.GetCol());
	for(unsigned i=0;i<right.GetRow();i++)
		for(unsigned j=0;j<right.GetCol();j++)
			temp(i,j)=right(i,j)*coef;
	return temp;
}

Matrix Matrix::transpose(void)const{
	Matrix temp(rep->col,rep->row);
	for(unsigned i=0;i<rep->row;i++)
		for(unsigned j=0;j<rep->col;j++)
			temp(j,i)=(*this)(i,j);
	return temp;
}

double Matrix33::I1(void)const{
	return (*this)(0,0)+(*this)(1,1)+(*this)(2,2);
}

double Matrix33::I2(void)const{
	double sum=0;
	for(unsigned i=0;i<3;i++)
	for(unsigned j=0;j<3;j++)
		sum+=(*this)(i,j)*(*this)(j,i);
	double I1=this->I1();
	return (I1*I1-sum)*.5;
}

double Matrix33::I3(void)const{
	return
		(*this)(0,0)*(*this)(1,1)*(*this)(2,2) +
		(*this)(0,1)*(*this)(1,2)*(*this)(2,0) +
		(*this)(0,2)*(*this)(1,0)*(*this)(2,1) -
		(*this)(0,2)*(*this)(1,1)*(*this)(2,0) -
		(*this)(0,0)*(*this)(1,2)*(*this)(2,1) -
		(*this)(0,1)*(*this)(1,0)*(*this)(2,2);
}

void Matrix33::identity(void){
	zero_out();
	for(unsigned i=0;i<3;i++)
		(*this)(i,i)=1;
}

void Matrix33::dev(void){
	double trace3=I1()/3.;
	for(unsigned i=0;i<3;i++)
		(*this)(i,i)-=trace3;
}

Matrix33 Matrix33::inv(void)const{
	double mdet=I3();
	Matrix temp(GetRow(),GetCol());
	temp(0,0)=((*this)(1,1)*(*this)(2,2) - (*this)(1,2)*(*this)(2,1))/mdet;
	temp(0,1)=((*this)(0,2)*(*this)(2,1) - (*this)(0,1)*(*this)(2,2))/mdet;
	temp(0,2)=((*this)(0,1)*(*this)(1,2) - (*this)(0,2)*(*this)(1,1))/mdet;
	temp(1,0)=((*this)(1,2)*(*this)(2,0) - (*this)(1,0)*(*this)(2,2))/mdet;
	temp(1,1)=((*this)(0,0)*(*this)(2,2) - (*this)(0,2)*(*this)(2,0))/mdet;
	temp(1,2)=((*this)(0,2)*(*this)(1,0) - (*this)(0,0)*(*this)(1,2))/mdet;
	temp(2,0)=((*this)(1,0)*(*this)(2,1) - (*this)(1,1)*(*this)(2,0))/mdet;
	temp(2,1)=((*this)(0,1)*(*this)(2,0) - (*this)(0,0)*(*this)(2,1))/mdet;
	temp(2,2)=((*this)(0,0)*(*this)(1,1) - (*this)(0,1)*(*this)(1,0))/mdet;
	return temp;
}

/*
int main(void){
	Matrix A(2,2);
	A(0,0)=3.1;
	A(0,1)=2.1;
	A(1,0)=3.4;
	A(1,1)=3.0;

	Matrix B(2,2);
	B(0,0)=1.1;
	B(0,1)=1.1;
	B(1,0)=1.4;
	B(1,1)=1.0;

	Matrix C(2,2);
	C(0,0)=1.1;
	C(0,1)=1.1;
	C(1,0)=1.4;
	C(1,1)=1.0;

	Matrix D;

	cout << "A"<<endl;
	A.print(cout);
	cout << "B"<<endl;
	B.print(cout);
	A=A+2.0*B*2.0;
	cout << "A"<<endl;
	A.print(cout);
	//D=A+B+C;
	//D.print(cout);

	cout << "33" << endl;
	Matrix33 X,Y,Z;
	X=Y+Z;
	return 0;
}
*/


#undef _DEBUG
