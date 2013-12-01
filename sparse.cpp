/*
ver.03.03
file sparse.cpp
*/

#include "sparse.h"
using std::vector;

const double Sparse_wiz::tolerance=1.0E-14;

Sparse_wiz::Sparse_wiz(bool symmetric):
N(0),
ia(NULL),
ja(NULL),
iap1(NULL),
jap1(NULL),
values(NULL),
isInitialized(false),
vja(NULL),
Sparse(symmetric)
{}

Sparse_wiz::~Sparse_wiz(void){
	if(values)
		delete[] values;
	if(ia)
		delete[] ia;
	if(ja)
		delete[] ja;
	if(vja)
		delete[] vja;
	if(iap1)
		delete[] iap1;
	if(jap1)
		delete[] jap1;
}

void Sparse_wiz::initialize(int N){
	if(isInitialized)
		return;
	isInitialized=true;
	this->N=N;
	ia=new int[N+1];
	vja=new vector<int>[N];
	for(int i=0;i<N;i++){
		ia[i]=i;
		vja[i].push_back(i);
	}
	ia[N]=N;
}

#ifdef __INTEL_COMPILER
void Sparse_wiz::solve_mkl(double* rhs,double* solution){
	if(!symmetric){
		for(unsigned row=0;row<N;row++)
			for(unsigned j=ia[row];j<ia[row+1];j++){
				unsigned col=ja[j];
				if(row<col)
					values[getSparseIndex(col,row)]=values[j];
			}
	}
	pardiso(rhs,solution,N,iap1,jap1,values);
}
#endif

void Sparse_wiz::add(int row,int col){
	if(row==col){
		return;
	}else if(vja[row].front()>col){
		//vja[row].push_front(col);
		vja[row].insert(vja[row].begin(),col);
	}else{
		vector<int>::iterator pja=vja[row].end();
		for(;pja!=vja[row].begin();){
			pja--;
			if(*pja==col)
				return;
			if(*pja<col)
				break;
		}
		vja[row].insert(++pja,col);
	}
	//fixme: this can be made faster by using frequency occurance
	for(int i=row+1;i<N+1;i++)
		ia[i]++;

	if(!symmetric && row<col)
		add(col,row);

}

void Sparse_wiz::add(int row,int col,double value){
	SPARSEINDEX si=getSparseIndex(row,col);
#ifdef _DEBUG
	if(si==-1){
		cout << "error in Sparse_wiz::add(value)" << endl;
		cout << row << " " << col << endl;
		return;
	}
#endif
	values[si]+=value;
}

double& Sparse_wiz::operator()(int row,int col){
	SPARSEINDEX si=getSparseIndex(row,col);
	return values[si];
}

void Sparse_wiz::allocate(void){
	vja_size=0;
	for(int i=0;i<N;i++)
		vja_size+=vja[i].size();

	ja=new int[vja_size];
	values=new double[vja_size];
	ia[N]=vja_size;

	int cnt=0;
	for(int i=0;i<N;i++)
		for(vector<int>::iterator pja=vja[i].begin();pja!=vja[i].end();pja++)
			ja[cnt++]=*pja;

	for(unsigned int i=0;i<vja_size;i++)
		values[i]=0.0;

	iap1=new int[N+1];
	for(int i=0;i<N+1;i++)
		iap1[i]=ia[i]+1;
	jap1=new int[vja_size];
	for(int i=0;i<vja_size;i++)
		jap1[i]=ja[i]+1;
}

void Sparse_wiz::zeroOut(void){
	for(unsigned int i=0;i<vja_size;i++)
		values[i]=0.0;
}

void Sparse_wiz::print(const char filename[]){
	using std::ofstream;
	using std::endl;

	ofstream fout(filename);

	for(MATRIXINDEX i=0;i<N;i++){
	for(MATRIXINDEX j=0;j<N;j++){
		SPARSEINDEX si=getSparseIndex(i,j);
		if(si==-1){
			fout << 0 << ",";
		}else{
			fout << values[si] << ",";
		}
	}
	fout << endl;
	}
}

int Sparse_wiz::getSparseIndex(MATRIXINDEX row,MATRIXINDEX col){
	int cnt=ia[row];
	for(vector<int>::iterator pja=vja[row].begin();pja!=vja[row].end();pja++){
		if(*pja==col)
			return cnt;
		if(*pja>col)
			return -1;
		cnt++;
	}
	return -1;
}


