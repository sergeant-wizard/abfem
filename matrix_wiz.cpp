/*
ver.03.03
file matrix_wiz.cpp
*/

#include "matrix_wiz.h"
using std::ofstream;
using std::endl;
using std::cout;
using std::string;
using std::ostream;
using std::left;
using std::setw;

const double MATRIX_WIZ::m_eps=1.0E-12;

MATRIX_WIZ::MATRIX_WIZ(void):
row(1),
col(1),
size(1),
data(NULL),
data_cnt(NULL),
input_flag(false)
{
//	cout << "constructor called" << endl;
}

MATRIX_WIZ::MATRIX_WIZ(unsigned _row,unsigned _col):
row(_row),
col(_col),
size(_row*_col),
input_flag(false)
{
//	cout << "constructor called" << endl;
	data=new double[size];
	if(!data)
		cout<<"allocation error in MATRIX_WIZ::MATRIX_WIZ()"<<endl;

	for(unsigned i=0;i<size;i++)
		data[i]=0.0;
}

MATRIX_WIZ::MATRIX_WIZ(MATRIX_WIZ& src):
row(src.row),
col(src.col),
size(src.size),
input_flag(false)
{
	data=new double[src.size];
	if(!data)
		cout<<"allocation failure in MATIRX::MATRIX_WIZ(MATRIX_WIZ&)"<<endl;
	for(unsigned i=0;i<size;i++)
		data[i]=src.data[i];
}

MATRIX_WIZ::MATRIX_WIZ(MATRIX_WIZ* src):
row(src->row),
col(src->col),
size(src->size),
input_flag(false)
{
	data=new double[src->size];
	if(!data)
		cout<<"allocation failure in MATIRX::MATRIX_WIZ(MATRIX_WIZ&)"<<endl;
	for(unsigned i=0;i<size;i++)
		data[i]=src->data[i];
}

MATRIX_WIZ::~MATRIX_WIZ(void){
	if(data)
		delete[] data;
	data=NULL;
}

void MATRIX_WIZ::Renew(unsigned row,unsigned col){
	this->row=row;
	this->col=col;
	size=row*col;

	if(data)
		delete[] data;
	data=new double[size];
	if(!data)
		cout << "allocation error in MATRIX_WIZ::Renew()" << endl;
	ZeroOut();
}

void MATRIX_WIZ::ZeroOut(void){
	if(!data)return;
	for(unsigned i=0;i<size;i++)
		data[i]=0.0;
}

bool MATRIX_WIZ::IsSqr(void){
	if(row==col)
		return true;
	else
		return false;
}

void MATRIX_WIZ::Input(unsigned num,...){
	va_list list;
	va_start(list,num);
	for(unsigned i=0;i<num;i++)
		data[i]=va_arg(list,double);
	va_end(list);
}

void MATRIX_WIZ::Divide(MATRIX_WIZ& smaller,unsigned top,unsigned bottom,unsigned left,unsigned right){
	unsigned srow=bottom-top+1;
	unsigned scol=right-left+1;
	for(unsigned i=0;i<srow;i++)
		for(unsigned j=0;j<scol;j++)
			smaller(i,j)=(*this)(top+i,left+j);
}

void MATRIX_WIZ::Identity(void){
#ifdef _DEBUG
	if(!IsSqr()){
		cout << "non-square matrix given to MATRIX_WIZ::Identity()" << endl;
		return;
	}
#endif
	unsigned cnt=0;
	for(unsigned i=0;i<row;i++){
		for(unsigned j=0;j<col;j++){
			if(i==j)
				data[cnt++]=1.0;
			else
				data[cnt++]=0.0;
		}
	}
}

void MATRIX_WIZ::Print(string file){
	if(!data)return;

	ofstream fout(file.c_str());
	unsigned cnt=0;
	for(unsigned i=0;i<row;i++){
		for(unsigned j=0;j<col;j++){
			fout << data[cnt++] << ",";
		}
		fout << endl;
	}
}

void MATRIX_WIZ::Eigen(MATRIX_WIZ& eigenvectors,MATRIX_WIZ& eigenvalues){
#ifdef _DEBUG
	if(!IsSqr()){
		cout << "error in Matrix::Eigen()" << endl;
		return;
	}
	eigenvectors.Renew(row,col);
	eigenvalues.Renew(row,1);
#endif

	unsigned iteration_cnt=1000;
	MATRIX_WIZ R(row,col);
	MATRIX_WIZ A(*this);
	eigenvectors.Identity();

	while(1){
		if(!(iteration_cnt--)){
			cout << "failed to calculate Matrix::Eigen()" << endl;
			return;
		}
		double maxa_ij=fabs(A(0,1));
		unsigned max_i=0,max_j=1;
		for(unsigned i=0;i<row;i++){
			for(unsigned j=i+1;j<col;j++){
				if(fabs(A(i,j))>maxa_ij){
					maxa_ij=fabs(A(i,j));
					max_i=i;
					max_j=j;
				}
			}
		}

		if(fabs(maxa_ij)<m_eps)
			break;

		double a_ii=A(max_i,max_i);
		double a_jj=A(max_j,max_j);
		double a_ij=A(max_i,max_j);
		double theta=atan2(2.0*a_ij,a_ii-a_jj)/2.0;

		R.Identity();
		R(max_i,max_i)= cos(theta);
		R(max_i,max_j)=-sin(theta);
		R(max_j,max_i)= sin(theta);
		R(max_j,max_j)= cos(theta);

		A=R.Transpose()*A*R;
		eigenvectors=eigenvectors*R;
	}
	for(unsigned i=0;i<row;i++)
		eigenvalues(i,0)=A(i,i);
}

double MATRIX_WIZ::Trace(void){
#ifdef _DEBUG
	if(!IsSqr()){
		cout << "error in Matrix::Trace()" << endl;
		return 0.0;
	}
#endif
	double sum=0;
	for(unsigned i=0;i<row;i++)
		sum+=data[index(i,i)];
	return sum;
}

double MATRIX_WIZ::Norm2(void){
	double sum=0.0;
	for(double* pdata=data;pdata-data<size;pdata++)
		sum+=(*pdata)*(*pdata);
	return sum;
}

void MATRIX_WIZ::Inv(MATRIX_WIZ& dest){
#ifdef _DEBUG
	if(!IsSqr()){
		cout << "error: non-square matrix given for Matrix::inv()" << endl;
		return;
	}
#endif
	MATRIX_WIZ buffer(*this);
	dest.Renew(row,col);
	dest.Identity();
	double diagonal;
	unsigned i,j,k;
	double temp;

	//forward
	for(i=0;i<row;i++){
		diagonal=buffer(i,i);
		if(diagonal==0.0){
			for(j=i+1;j<row;j++){
				if(buffer(j,i)!=0.0){
					buffer.gauss(i,j);
					dest.gauss(i,j);
					diagonal=buffer(i,i);
					break;
				}
			}
		}
		buffer.gauss(i,1.0/diagonal);
		dest.gauss(i,1.0/diagonal);
		for(j=i+1;j<row;j++){
			if(temp=-buffer(j,i)){
				for(k=i;k<col;k++){
					buffer.data[index(j,k)]+=buffer.data[index(i,k)]*temp;
//					if(fabs(buffer.data[index(j,k)])<m_eps)
//						buffer.data[index(j,k)]=0.0;
				}
				dest.gauss(j,i,temp);
			}
		}
	}

	//backward
	for(i=row;i-->1;){
		for(j=i;j-->0;){
			if(temp=-buffer.data[index(j,i)]){
				for(k=i;k<col;k++){
					buffer.data[index(j,k)]+=buffer.data[index(i,k)]*temp;
//					if(fabs(buffer.data[index(j,k)])<m_eps)
//						buffer.data[index(j,k)]=0.0;
				}
				dest.gauss(j,i,temp);
			}
		}
	}
	return;
}

MATRIX_WIZ& MATRIX_WIZ::Transpose(void){
	static MATRIX_WIZ ret;
	ret.Renew(col,row);
	double* psrc=data;
	double* pdest=ret.data;

	for(unsigned i=0;i<col;i++){
		for(unsigned j=0;j<row;j++){
			*pdest=(*this)(j,i);
			pdest++;
		}
	}
	return ret;
}

MATRIX& MATRIX_WIZ::operator=(MATRIX_WIZ& right){
#ifdef _DEBUG
	if(row!=right.row || col!=right.col){
		cout << "invalid operation for MATRIX_WIZ::=" << endl;
		return *this;
	}
	if(!data || !right.data){
		cout << "invalid matrix for MATRIX_WIZ::=" << endl;
		return *this;
	}
#endif
	for(unsigned i=0;i<size;i++)
		data[i]=right[i];
	return *this;
}

MATRIX& MATRIX_WIZ::operator+=(MATRIX& right){
#ifdef _DEBUG
	if(row!=right.row || col!=right.col){
		cout << "invalid operation for MATRIX_WIZ::+=" << endl;
		return *this;
	}
	if(!data || !right.data){
		cout << "invalid matrix for MATRIX_WIZ::+=" << endl;
		return *this;
	}
#endif

	for(unsigned i=0;i<size;i++)
		data[i]+=right[i];
	return *this;
}

MATRIX& MATRIX_WIZ::operator-=(MATRIX& right){
#ifdef _DEBUG
	if(row!=right.row || col!=right.col){
		cout << "invalid operation for MATRIX_WIZ::-=" << endl;
		return *this;
	}
	if(!data || !right.data){
		cout << "invalid matrix for MATRIX_WIZ::-=" << endl;
		return *this;
	}
#endif

	for(unsigned i=0;i<size;i++)
		data[i]-=right[i];
	return *this;
}

MATRIX& MATRIX_WIZ::operator*=(double coef){
	if(!data)
		return *this;
	for(unsigned i=0;i<size;i++)
		data[i]*=coef;
	return *this;
}

MATRIX& MATRIX_WIZ::operator*=(MATRIX& right){
	if(!data){
		cout << "error in MATRIX_WIZ::operator*=(MATRIX_WIZ)" << endl;
		return *this;
	}
	MATRIX_WIZ temp=*this;
	for(unsigned i=0;i<row;i++)
	for(unsigned j=0;j<col;j++){
		double sum=0;
		for(unsigned k=0;k<row;k++)
			sum+=temp(i,k)*right(k,j);
		(*this)(i,j)=sum;
	}
	return *this;
}

ostream& operator<<(ostream& os,MATRIX_WIZ& matrix){
	if(!matrix.data)
		return os;
	cout << endl;
	unsigned cnt=0;
	for(unsigned i=0;i<matrix.row;i++){
		for(unsigned j=0;j<matrix.col;j++){
			cout << left << setw(MATRIX_WIZ::width) << matrix.data[cnt++] << "\t";
		}
		cout << endl;
	}
	return os;
}

MATRIX_WIZ& operator+(MATRIX_WIZ& left,MATRIX_WIZ& right){
	static MATRIX_WIZ ret(left.row,left.col);
#ifdef _DEBUG
	if(left.row!=right.row || left.col!=right.col){
		cout << "invalid operation for MATRIX_WIZ::+" << endl;
		return ret;
	}
	if(!left.data || !right.data){
		cout << "invalid matrix for MATRIX_WIZ::+" << endl;
		return ret;
	}
#endif
	for(unsigned i=0;i<left.size;i++)
		ret.data[i]=left.data[i]+right.data[i];
	return ret;
}

MATRIX_WIZ& operator-(MATRIX_WIZ& left,MATRIX_WIZ& right){
	static MATRIX_WIZ ret(left.row,left.col);
#ifdef _DEBUG
	if(left.row!=right.row || left.col!=right.col){
		cout << "invalid operation for MATRIX_WIZ::-" << endl;
		return ret;
	}
	if(!left.data || !right.data){
		cout << "invalid matrix for MATRIX_WIZ::-" << endl;
		return ret;
	}
#endif
	for(unsigned i=0;i<left.size;i++)
		ret.data[i]=left.data[i]-right.data[i];
	return ret;
}

MATRIX_WIZ& operator*(MATRIX_WIZ& mat1,MATRIX_WIZ& mat2){
	static MATRIX_WIZ ret;

#ifdef _DEBUG
	if(mat1.col!=mat2.row){
		cout << "error in MATRIX_WIZ::*()" << endl;
		return ret;
	}
#endif
	MATRIX_WIZ mat1_(mat1);
	MATRIX_WIZ mat2t(mat2.col,mat2.row);
	mat2t=mat2.Transpose();
	ret.Renew(mat1.row,mat2.col);
	unsigned cnt=0;

	for(unsigned i=0;i<mat1_.row;i++){
		for(unsigned j=0;j<mat2.col;j++){
			ret.data[cnt]=0.0;
			for(unsigned k=0;k<mat1_.col;k++)
				ret.data[cnt]+=mat1_(i,k)*mat2t(j,k);
			cnt++;
		}
	}
	return ret;
}

double& MATRIX_WIZ::operator()(unsigned i,unsigned j){
	if(!data)
		return data[0];
	else
		return data[i*col+j];
}

double MATRIX_WIZ::operator()(unsigned i,unsigned j)const{
	if(!data)
		return data[0];
	else
		return data[i*col+j];
}

double MATRIX_WIZ::Atm(unsigned rown,unsigned coln){
	if(!data)return 0.0;
	return data[(rown-1)*col+coln-1];
}

unsigned MATRIX_WIZ::index(unsigned _row,unsigned _col){
	return _row*col+_col;
}

void MATRIX_WIZ::gauss(unsigned _row,double coef){
	for(unsigned i=0;i<col;i++)
		if(data[index(_row,i)])
			data[index(_row,i)]*=coef;
}

void MATRIX_WIZ::gauss(unsigned row1,unsigned row2,double coef){
	for(unsigned i=0;i<col;i++)
		if(data[index(row2,i)])
			data[index(row1,i)]+=data[index(row2,i)]*coef;
}

void MATRIX_WIZ::gauss(unsigned row1,unsigned row2){
	unsigned i;
	double buffer;
	for(i=0;i<col;i++){
		buffer=data[index(row1,i)];
		data[index(row1,i)]=data[index(row2,i)];
		data[index(row2,i)]=buffer;
	}
}

//MATRIX33_WIZ
MATRIX33_WIZ::MATRIX33_WIZ(){
	Renew(3,3);
}

MATRIX33_WIZ::MATRIX33_WIZ(MATRIX33& src){
	Renew(3,3);
#ifdef _DEBUG
	if(src.row!=3 || src.col!=3){
		cout << "error in MATRIX33_WIZ::MATRIX33_WIZ" << endl;
		return;
	}
#endif
	for(unsigned i=0;i<9;i++)
		data[i]=src[i];
		
}

MATRIX33_WIZ::MATRIX33_WIZ(MATRIX33* src){
	Renew(3,3);
#ifdef _DEBUG
	if(src->row!=3 || src->col!=3){
		cout << "error in MATRIX33_WIZ::MATRIX33_WIZ" << endl;
		return;
	}
#endif
	for(unsigned i=0;i<9;i++)
		data[i]=(*src)[i];
}

void MATRIX33_WIZ::Inv(MATRIX33& dest){
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
}

MATRIX& MATRIX33_WIZ::operator=(MATRIX33_WIZ& right){
#ifdef _DEBUG
	if(right.col!=3 || right.col!=3){
		cout << "invalid operation for MATRIX33_WIZ::=" << endl;
		return *this;
	}
#endif
	for(unsigned i=0;i<9;i++)
		data[i]=right[i];
	return *this;
}

MATRIX& MATRIX33_WIZ::operator=(MATRIX_WIZ& right){
#ifdef _DEBUG
	if(right.col!=3 || right.col!=3){
		cout << "invalid operation for MATRIX33_WIZ::=" << endl;
		return *this;
	}
#endif
	for(unsigned i=0;i<9;i++)
		data[i]=right[i];
	return *this;
}

MATRIX33_WIZ& MATRIX33_WIZ::Transpose(void){
	static MATRIX33_WIZ ret;
	ret.Input(9,
		(*this)(0,0),
		(*this)(1,0),
		(*this)(2,0),
		(*this)(0,1),
		(*this)(1,1),
		(*this)(2,1),
		(*this)(0,2),
		(*this)(1,2),
		(*this)(2,2)
	);
	return ret;
}

void MATRIX33_WIZ::dev(void){
	double trace3=I1()/3.;
	for(unsigned i=0;i<3;i++)
		(*this)(i,i)-=trace3;
}

/*obsolete
double MATRIX33_WIZ::det(void)const{
	return
		(*this)(0,0)*(*this)(1,1)*(*this)(2,2) +
		(*this)(0,1)*(*this)(1,2)*(*this)(2,0) +
		(*this)(0,2)*(*this)(1,0)*(*this)(2,1) -
		(*this)(0,2)*(*this)(1,1)*(*this)(2,0) -
		(*this)(0,0)*(*this)(1,2)*(*this)(2,1) -
		(*this)(0,1)*(*this)(1,0)*(*this)(2,2);
}

double MATRIX33_WIZ::I1(void)const{
	return (*this)(0,0)+(*this)(1,1)+(*this)(2,2);
}

double MATRIX33_WIZ::I2(void)const{
	MATRIX33_WIZ temp=*this;
	temp*=*this;
	double I1=this->I1();
	return (I1*I1-temp.I1())*0.5;
}
*/
