/*
ver.03.03
file element.h
*/

#pragma once
#include <vector>
#include "sparse.h"
#include "matrix.h"
#include "node.h"

class Element{
public:
	Element(void){};
	virtual ~Element(){};
	virtual void initialize(void)=0;
	virtual void initialize_sparse(Sparse& K)=0;
	virtual void update_loadstep(void)=0;
	virtual void update_iteration(void)=0;
	virtual void calcK(Sparse& K)=0;
	virtual void calcT(MATRIX& T)=0;
	virtual Node*& get_node(unsigned)=0;
	virtual void sync_Kmatrix(Sparse& K,bool symmetric=true)=0;

	template<class Parent,class Child>
	static void Up_cast(std::vector<Parent*>& parent,std::vector<Child*>& child){
		for(unsigned i=0;i<child.size();i++)
			parent.push_back(child[i]);
	};
};

class Element_output{
public:
	virtual ~Element_output(){};
	virtual double getStress(unsigned,unsigned,unsigned)const=0;
	virtual double get_F_(unsigned,unsigned,unsigned)const=0;
	virtual double get_J_(unsigned)const=0;
	virtual double get_pressure(unsigned)const=0;
	virtual double get_potential(unsigned)const=0;
	virtual Node*& get_node(unsigned)=0;

	virtual double get_volume(void)const{return 0;};

	//issue11
	virtual unsigned get_effective_flag(void)const{return 1;};
};

class Element_UP:virtual public Element{
public:
	Element_UP(const unsigned num_dof_p):
		num_dof_p(num_dof_p)
	{
		m_p=new double[num_dof_p];
	};

	virtual ~Element_UP(){
		delete[] m_p;
	};

	virtual void calcK_UP(Sparse& K)=0;
	virtual void calcT_p(MATRIX& T)=0;

	virtual void put_pressureID(unsigned pi)const{m_pressure_id=pi;};
	virtual void put_pressureID_(unsigned pi)const{m_pressure_id_=pi;};
	virtual unsigned get_pressureID(void)const{return m_pressure_id;};
	virtual unsigned get_pressureID_(void)const{return m_pressure_id_;};

	virtual void increment_pressure(unsigned alpha,double dp){
		m_p[alpha]+=dp;
	};
	//fixme
	virtual double getStress_hat(unsigned alpha,unsigned,unsigned)const{return 0;};

	//issue12
	virtual void calcBodyForce(MATRIX& T){};

protected:
	const unsigned num_dof_p;
	double* m_p;
	mutable unsigned m_pressure_id;
	mutable unsigned m_pressure_id_;
};

