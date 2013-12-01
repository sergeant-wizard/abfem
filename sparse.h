/*
ver.03.03
file sparse.h
*/

/*
+Sparse
|-Sparse_pardiso
|-Sparse_lapack
*/

#pragma once
#ifdef __INTEL_COMPILER
#include "mkl_pardiso.h"
#include "mkl_types.h"
#include <omp.h>
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <list>

class Sparse{
public:
	Sparse(void):ia(NULL),ja(NULL){};
	//Sparse(const int N):ia(NULL),ja(NULL){};
	Sparse(bool symmetric):symmetric(symmetric),ia(NULL),ja(NULL){};
	virtual ~Sparse(){
		if(ia)
			delete[] ia;
		if(ja)
			delete[] ja;
	};
	virtual void initialize(const int N)=0;
	virtual void add(int row,int col){};
	virtual void add(int row,int col,double value)=0;
	virtual double& operator()(int row,int col){
		std::cout<<"override me"<<std::endl;return dummy;
	};
	virtual void allocate(void)=0;
	virtual void zeroOut(void)=0;
	virtual void print(const char filename[])=0;
	virtual void solve(double* rhs,double* solution)=0;
protected:
	void pardiso(double* rhs,double* solution,int N,int* ia_,int* ja_,double* values,bool subtracted=false){
		using std::cout;
		using std::endl;
#ifdef __INTEL_COMPILER
		MKL_INT mtype;

		if(symmetric)
			mtype = -2; // Real symmetric matrix
		else
			mtype = 11; // Real non-symmetric matrix

		//MKL_INT mtype = 2; // Real symmetric P.D matrix

		if(!ia)
			ia=new int[N+1];
		if(!ja)
			ja=new int[ia_[N]];

		if(subtracted){
			for(int i=0;i<N+1;i++)
				ia[i]=ia_[i]+1;
			for(int i=0;i<ia_[N];i++)
				ja[i]=ja_[i]+1;
		}else{
			for(int i=0;i<N+1;i++)
				ia[i]=ia_[i];
			for(int i=0;i<ia_[N];i++)
				ja[i]=ja_[i];
		}

		MKL_INT n = N;
		MKL_INT nrhs = 1;  // Number of right hand sides
		void *pt[64]={NULL}; //Internal solver memory pointer pt
		//Pardiso control parameters
		MKL_INT iparm[64]={0};
		MKL_INT maxfct, mnum, phase, error, msglvl;
		// Auxiliary variables
		double ddum;  // Double dummy 
		MKL_INT idum; // Integer dummy
		/////////////////////////////////////////////////////////////////////////
		// ..  Setup Pardiso control parameters. 
		/////////////////////////////////////////////////////////////////////////
		iparm[0] = 1;  // No solver default
		iparm[1] = 2;  // Fill-in reordering from METIS
		//iparm[2] = 1;  // Numbers of processors, value of OMP_NUM_THREADS
		iparm[2] = omp_get_num_threads();
		iparm[3] = 0;  // No iterative-direct algorithm 
		iparm[4] = 0;  // No user fill-in reducing permutation */
		iparm[5] = 0;  // Write solution into x */
		iparm[6] = 0;  // Not in use */
		iparm[7] = 2;  // Max numbers of iterative refinement steps */
		iparm[8] = 0;  // Not in use */
		iparm[9] = 13; // Perturb the pivot elements with 1E-13 */
		iparm[10] = 1; // Use nonsymmetric permutation and scaling MPS */
		iparm[11] = 0; // Not in use */
		iparm[12] = 0; // Maximum weighted matching algorithm is switched-off 
					   // (default for symmetric). 
					   // Try iparm[12] = 1 in case of inappropriate accuracy 
		iparm[13] = 0; // Output: Number of perturbed pivots */
		iparm[14] = 0; // Not in use */
		iparm[15] = 0; // Not in use */
		iparm[16] = 0; // Not in use */
		iparm[17] = -1;// Output: Number of nonzeros in the factor LU
		iparm[18] = -1;// Output: Mflops for LU factorization
		iparm[19] = 0; // Output: Numbers of CG Iterations */
		maxfct = 1;    // Maximum number of numerical factorizations
		mnum = 1;      // Which factorization to use.
//#ifdef _DEBUG
		//msglvl = 1;    // Print statistical information in file
//#else
		msglvl = 0;    // n'est-pa print
//#endif
		error = 0;     // Initialize error flag

		/////////////////////////////////////////////////////////////////////////
		// .. Reordering and Symbolic Factorization. This step also allocates
		// all memory that is necessary for the factorization.
		/////////////////////////////////////////////////////////////////////////
		phase = 11;
		PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
			&n, values, ia, ja, &idum, &nrhs,
			iparm, &msglvl, &ddum, &ddum, &error);
		if (error) {
			cout << "ERROR during symbolic factorization: " << error << endl;
			exit(1);
		}
		/////////////////////////////////////////////////////////////////////////
		// .. Numerical factorization.
		/////////////////////////////////////////////////////////////////////////
		phase = 22;
		PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
			&n, values, ia, ja, &idum, &nrhs,
			iparm, &msglvl, &ddum, &ddum, &error);
		if (error) {
			cout << "ERROR during numerical factorization: " << error << endl;
			exit(2);
		}
		/////////////////////////////////////////////////////////////////////////
		// .. Back substitution and iterative refinement.
		/////////////////////////////////////////////////////////////////////////
		phase = 33;
		iparm[7] = 2;  //Max numbers of iterative refinement steps. 
		PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
			&n, values, ia, ja, &idum, &nrhs,
			iparm, &msglvl, rhs, solution, &error);
		if (error) {
			cout << "ERROR during solution: " << error << endl;
			exit(3);
		}
		/////////////////////////////////////////////////////////////////////////
		// .. Termination and release of memory.
		/////////////////////////////////////////////////////////////////////////
		phase = -1; // Release internal memory
		PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
			&n, &ddum, ia, ja, &idum, &nrhs,
			iparm, &msglvl, &ddum, &ddum, &error);
#endif
	}
protected:
	bool symmetric;
private:
	int* ia;
	int* ja;
	double dummy;
};

class Sparse_wiz:public Sparse{
	typedef int FULLINDEX;
	typedef int SPARSEINDEX;
	typedef int MATRIXINDEX;
public:
	Sparse_wiz(bool symmetric=false);
	~Sparse_wiz();
	void initialize(const int N);
	void solve(double* rhs,double* solution){
//issue6
#ifdef __INTEL_COMPILER
		solve_mkl(rhs,solution);
#else
		;
#endif
	};
	void add(int row,int col);
	void add(int row,int col,double value);
	double& operator()(int row,int col);
	void allocate(void);
	void zeroOut(void);
	void print(const char filename[]);
	int size(void)const{
		return N;
	};
private:
	void solve_mkl(double* rhs,double* solution);
	int getSparseIndex(MATRIXINDEX row,MATRIXINDEX col);

	const static double tolerance;
	int N;
	std::vector<int> *vja;
	int *iap1,*jap1;
	int *ia,*ja;
	double *values;
	bool isInitialized;
	int vja_size;
};


