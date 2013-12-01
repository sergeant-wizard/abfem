/*
ver.03.03
file material.h
*/

/*
0 bonet neo bookean imcompressible
1 bonet neo hookean compressible
2 owen neo hookean
3 St. Venant
4 simo neo-hookean
*/

#pragma once
#include "common.h"

template <const unsigned DIM,class F_TYPE=void*>
class Elastic{
public:
	Elastic(const double E,const double nu):
		E(E),
		nu(nu),
		lambda(E*nu/(1.0+nu)/(1.0-2.0*nu)),
		mu(E*0.5/(1.0+nu)),
		kappa(E/3.0/(1.0-2.0*nu))
	{};
	
	virtual void getC_e(double dest[][DIM][DIM][DIM],F_TYPE* F)=0;
	virtual void getSigma(F_TYPE* dest,F_TYPE* F)=0;
	virtual double potential(const F_TYPE* F){return 0;};
protected:
    const double E;
	const double nu;
	const double lambda;
	const double mu;
	const double kappa;
};

//truly incompressible to be precise, otherwise...
template <const unsigned DIM,class F_TYPE=void*>
class Incompressible: virtual public Elastic<DIM,F_TYPE>{
public:
	Incompressible(const double E):Elastic<DIM,F_TYPE>(E,0.5){};
};

template <const unsigned DIM,class F_TYPE=void*>
class Plastic: public Elastic<DIM,F_TYPE>{
public:
	Plastic(const double E,const double nu,const double tau0=0.0,const double H=0.0):
	Elastic<DIM,F_TYPE>(E,nu),
	tau0(tau0),
	H(H)
	{};
protected:
	const double tau0;
	const double H;
};

