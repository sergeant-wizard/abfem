/*
ver.03.03
file common.h
*/

#pragma once
#include <iostream>
#include <string>
#include <fstream>

#ifndef _MSC_VER
#define	sprintf_s(buffer,buffer_size,stringbuffer,...) (sprintf(buffer,stringbuffer,__VA_ARGS__))
#define	sscanf_s(src,format,...) (sscanf(src,format,__VA_ARGS__))
#endif

//#define _DEBUG
#define CHECK_MARKER
#define SHORT
//#define ALPHA
//#define MAX_VOLUME
#define ITERATOR
#define FAST_ISINSIDE
//#define TOTAL_LAGRANGE
//#define MASS

//issue11
//#define BETA

namespace ABFEM{
namespace CONSTANTS{
const double PI= 3.1415926535897932;
const double PI2=6.2831853071795864;
const double DIM=3;
static const double g_eps=1.0E-13;
static const double g_geo_eps=1.0E-9;
static const double g_diffuse=1.0E-3;
static const int g_VP_div=10;
static const double g_kpp_diag=1E-10;
static const double g_mass_coef=0;
//static const double g_mass_coef=1.0;
};

namespace MATH{
template <class Type> Type SQR(Type arg){
	return arg*arg;
}
template <class Type> inline const Type& MAX(const Type& a,const Type& b){
	return (a>b)?a:b;
}
template <class Type> inline const Type& MIN(const Type& a,const Type& b){
	return (a<b)?a:b;
}
const static double kronecker[3][3]={
	{1.0,0.0,0.0},
	{0.0,1.0,0.0},
	{0.0,0.0,1.0}
};

const static double kronecker2[3][3][3][3]={
	{//i=0
		{
			{1.0,0.0,0.0},
			{0.0,1.0,0.0},
			{0.0,0.0,1.0}
		},{0.0},{0.0}
	},{//i=1
		{0.0},{
			{1.0,0.0,0.0},
			{0.0,1.0,0.0},
			{0.0,0.0,1.0}
		},{0.0}
	},{//i=2
		{0.0},{0.0},{
			{1.0,0.0,0.0},
			{0.0,1.0,0.0},
			{0.0,0.0,1.0}
		}
	}
};
};

template <class Type> void SWAP(Type& a,Type& b){
	Type tmp=a;
	a=b;
	b=tmp;
}

template <class Type> void SORT(int num,Type* arg){
	for(int i=0;i<num-1;i++){
		Type min=arg[i];
		int min_index=i;
		for(int j=i+1;j<num;j++)
			if(arg[j]<min){
				min_index=j;
				min=arg[j];
			}
		SWAP(arg[i],arg[min_index]);
	}
}

template <class Type> void SAFE_NEW(Type*& arg,unsigned int size){
	if(arg)
		delete[] arg;
	arg=new Type[size];
}

template <class Type> void SAFE_DELETE(Type*& arg){
	if(arg)
		delete[] arg;
	arg=NULL;
}


const static double alternate[3][3][3]={
	{
		{0.0},
		{0.0,0.0,1.0},
		{0.0,-1.0,0.0}
	},{
		{0.0,0.0,-1.0},
		{0.0},
		{1.0,0.0,0.0}
	},{
		{0.0,1.0,0.0},
		{-1.0,0.0,0.0},
		{0.0}
	}
};


static std::string gstr;

class Progress{
public:
	Progress(int all):all(all),current(-1){};
	void display(int step){
		using std::cout;
		if(step==all-1){
			cout << "\r\x1b[K";
			cout << "[";
			for(int i=0;i<div;i++)
				cout << "*";
			cout << "]\n";
		}else if((step*div)/all==current){
			return;
		}else{
			current++;
			cout << "\r\x1b[K";
			cout << "[";
			for(int i=0;i<current;i++)
				cout << "*";
			for(int i=current;i<div;i++)
				cout << "-";
			cout << "]";
		}
		cout.flush();
	};
private:
	static const int div=10;
	const int all;
	int current;
};

#ifdef _DEBUG
#include "matrix_wiz.h"
static void printc(double m_c[][3][3][3]){
	using std::cout;
	MATRIX_WIZ D(6,6);
	D.Input(36,
		m_c[0][0][0][0],m_c[0][0][1][1],m_c[0][0][2][2],m_c[0][0][0][1],m_c[0][0][0][2],m_c[0][0][1][2],
		m_c[1][1][0][0],m_c[1][1][1][1],m_c[1][1][2][2],m_c[1][1][0][1],m_c[1][1][0][2],m_c[1][1][1][2],
		m_c[2][2][0][0],m_c[2][2][1][1],m_c[2][2][2][2],m_c[2][2][0][1],m_c[2][2][0][2],m_c[2][2][1][2],
		m_c[0][1][0][0],m_c[0][1][1][1],m_c[0][1][2][2],m_c[0][1][0][1],m_c[0][1][0][2],m_c[0][1][1][2],
		m_c[0][2][0][0],m_c[0][2][1][1],m_c[0][2][2][2],m_c[0][2][0][1],m_c[0][2][0][2],m_c[0][2][1][2],
		m_c[1][2][0][0],m_c[1][2][1][1],m_c[1][2][2][2],m_c[1][2][0][1],m_c[1][2][0][2],m_c[1][2][1][2]
	);
	cout<<D;
}
#endif

};

