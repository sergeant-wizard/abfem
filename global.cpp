/*
ver.03.03
file global.cpp
*/

#include "global.h"
using std::vector;
using std::cout;
using std::endl;
using std::setw;
using std::left;
using std::scientific;

//public
Global::Global(
	const unsigned loadSteps,
	const double tolerance_equivalence,
	const unsigned max_iteration,
	const unsigned interval,
	bool symmetric
):
loadSteps(loadSteps),
tolerance_equivalence(tolerance_equivalence),
max_iteration(max_iteration),
interval(interval),
symmetric(symmetric),
prepost(NULL),
K(NULL)
{
	K=new Default_sparse(symmetric);
}

Global::~Global(){
	if(prepost)
		delete prepost;
	prepost=NULL;
	if(K)
		delete K;
}

bool Global::solve(void){
	using namespace ABFEM::MATH;
	//issue12
	extern int g_loadstep_cnt;
	for(unsigned load_step=0;load_step<loadSteps;load_step++){
		//issue12
		g_loadstep_cnt=load_step+1;
		#ifdef SHORT
		cout << "\r\x1b[K";
		#endif
		cout << "load step:" << load_step << "/" << loadSteps << endl;
		int iteration_cnt=0;
		double tol;

		update_loadstep();
		update_iteration();

		int debug_cnt=0;
	
		do{
			#ifdef SHORT
			cout << "\x1b[K";
			#endif
			cout << "iteration:" << iteration_cnt << endl;
			iteration();

			tol=Rminus.Norm2()/(T.Norm2()+F.Norm2());

			#ifdef SHORT
			cout << "\x1b[K";
			#endif
			cout << "tolerance:" << left << scientific << setw(16) << sqrt(tol) << endl;
			if(
				isnan(tol) || isinf(tol) || tol<0.0 ||
				iteration_cnt++==max_iteration
			){
				cout << endl << "####didn't converge" << endl;
				if(((load_step)%(loadSteps/interval)))
					prepost->post_process(load_step-1);
				return false;
			}
			#ifdef SHORT
			cout << "\x1b[2A\r";
			#endif
			cout.flush();

   			//prepost->post_process(debug_cnt++);


		}while(tol>SQR(tolerance_equivalence));
		after_iteration(load_step);
		#ifdef SHORT
		cout << "\x1b[2A\r";
		#endif
		cout << "####" << iteration_cnt << " iterations" << endl;

	}
	#ifdef SHORT
	cout << "\x1b[4B";
	#endif
	return true;
}

//protected
void Global::after_preprocess(void){
	find_matrix_size();

	for(unsigned int i=0;i<velements.size();i++)
		velements[i]->initialize();

	//global vectors, matrices
	for(int i=0;i<num_dof_dirichlet+num_dof_nond;i++)
		vdof[i]->u=0.0;

	T.Renew(num_row+num_dof_dirichlet,1);
	F.Renew(num_row+num_dof_dirichlet,1);
	Rminus.Renew(num_row,1);

	K->initialize(num_row);
	ABFEM::Progress progress(velements.size());
	cout << "initializing sparse matrix..." << endl;
	for(unsigned int i=0;i<velements.size();i++){
		progress.display(i);
		velements[i]->initialize_sparse(*K);
	}
	cout << "...done" << endl << endl;
}

void Global::update_loadstep(void){
	for(vector<Dof*>::iterator it=vdof.begin();it!=vdof.end();it++){
		(*it)->X=(*it)->x;
		(*it)->u=0;
	}
	for(int i=0;i<num_dof_dirichlet;i++){
		vdof[i]->x+=vdof[i]->dirichlet/(double)loadSteps;
		vdof[i]->u=vdof[i]->dirichlet/(double)loadSteps;
	}
	for(unsigned i=0;i<velements.size();i++)
		velements[i]->update_loadstep();
}

void Global::update_iteration(void){
#ifdef __INTEL_COMPILER
	#pragma omp parallel for
#endif
	for(unsigned int i=0;i<velements.size();i++)
		velements[i]->update_iteration();
	findT();
}

void Global::after_iteration(unsigned load_step){
	prepost->save_result(load_step);
	if(!((load_step+1)%(loadSteps/interval)) || !load_step)
		prepost->post_process(load_step);
}

void Global::allocate_K(void){
	K->allocate();
	for(unsigned int i=0;i<velements.size();i++)
		velements[i]->sync_Kmatrix(*K);
}

void Global::findK(void){
	K->zeroOut();
	for(unsigned int i=0;i<velements.size();i++){
		velements[i]->calcK(*K);
	}
}

//private
void Global::iteration(void){
	findK();
	calculate();
	update_iteration();
}

void Global::findT(void){
	T.ZeroOut();
	for(unsigned int i=0;i<velements.size();i++){
		velements[i]->calcT(T);
	}
}

void Global::calculate(void){
	for(unsigned i=0;i<num_row;i++)
		Rminus(i,0)=F(i+num_dof_dirichlet,0)-T(i+num_dof_dirichlet,0);

	double* delta_u=new double[num_row];
	//issue6
#ifdef __MKL_PARDISO_H
	K->solve(Rminus.GetData(),delta_u);
#endif
	affect_solution(delta_u);
	delete[] delta_u;
}

//for(int i=0;i<num_comp_all;i++){
	//cout << left << setw(8) << vcinfo[i]->X << "\t" << vcinfo[i]->x << endl;
//}

