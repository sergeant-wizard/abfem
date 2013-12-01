/*
ver.03.03
file face.h
*/

#pragma once
#include <math.h>
#include "common.h"

class Face{
public:
	typedef double VECTOR[3];
	typedef VECTOR LINE[2];
	typedef VECTOR FACE[3];

	struct CutResult{
		double volume;
		double area;
	};

	Face(void):input_cnt(0){};
	Face(VECTOR normal,VECTOR point);
	void initialize(double src);
	double distance(VECTOR point);
	static bool Cross(Face& face1,Face& face);

	double getNormal(int component){return normal[component];};
	double* getVertice(int index){return (double*)&vertice[index];};

private:
	bool cross(double result[2],LINE line);
	VECTOR point;//=vertice[0]
	VECTOR normal;
	VECTOR vertice[3];
	int input_cnt;
	const static double m_eps;

	static double det(double src[9]){
		return
		src[0]*src[4]*src[8] +
		src[1]*src[5]*src[6] +
		src[2]*src[3]*src[7] -
		src[2]*src[4]*src[6] -
		src[0]*src[5]*src[7] -
		src[1]*src[3]*src[8];
	};

	static bool inv(double dest[9],double src[9]){
		double mdet=det(src);
		if(mdet<m_eps)
			return false;

		dest[0]=(src[4]*src[8] - src[5]*src[7])/mdet;
		dest[1]=(src[2]*src[7] - src[1]*src[8])/mdet;
		dest[2]=(src[1]*src[5] - src[2]*src[4])/mdet;
		dest[3]=(src[5]*src[6] - src[3]*src[8])/mdet;
		dest[4]=(src[0]*src[8] - src[2]*src[6])/mdet;
		dest[5]=(src[2]*src[3] - src[0]*src[5])/mdet;
		dest[6]=(src[3]*src[7] - src[4]*src[6])/mdet;
		dest[7]=(src[1]*src[6] - src[0]*src[7])/mdet;
		dest[8]=(src[0]*src[4] - src[1]*src[3])/mdet;
		return true;
	};

	static void getNormal(VECTOR dest,VECTOR v0,VECTOR v1,VECTOR v2){
		dest[0]=
			v0[1]*v1[2]-v0[2]*v1[1]+
			v1[1]*v2[2]-v1[2]*v2[1]+
			v2[1]*v0[2]-v2[2]*v0[1];
		dest[1]=
			v0[2]*v1[0]-v0[0]*v1[2]+
			v1[2]*v2[0]-v1[0]*v2[2]+
			v2[2]*v0[0]-v2[0]*v0[2];
		dest[2]=
			v0[0]*v1[1]-v0[1]*v1[0]+
			v1[0]*v2[1]-v1[1]*v2[0]+
			v2[0]*v0[1]-v2[1]*v0[0];
		normalize(dest);
	};

	static void outerProduct(VECTOR dest,VECTOR v0,VECTOR v1){
		dest[0]=v0[1]*v1[2]-v0[2]*v1[1];
		dest[1]=v0[2]*v1[0]-v0[0]*v1[2];
		dest[2]=v0[0]*v1[1]-v0[1]*v1[0];
	};

	static void normalize(VECTOR vec){
		using namespace ABFEM::MATH;
		double scale=SQR(vec[0])+SQR(vec[1])+SQR(vec[2]);
		for(int i=0;i<3;i++)
			vec[i]/=sqrt(scale);
	};

	static bool getTransform(double dest[9],VECTOR v0,VECTOR v1,VECTOR v2){
		double transform_inv[9]={
			v0[0],v1[0],v2[0],
			v0[1],v1[1],v2[1],
			v0[2],v1[2],v2[2]
		};
		return inv(dest,transform_inv);
	};


	static void getNormalTransform(double dest[9],VECTOR v0,VECTOR v1,VECTOR v2){
		VECTOR base0={v1[0]-v0[0],v1[1]-v0[1],v1[2]-v0[2]};
		normalize(base0);

		VECTOR base1;
		getNormal(base1,v0,v1,v2);

		VECTOR base2;
		outerProduct(base2,base0,base1);

		dest[0]=base0[0];
		dest[1]=base0[1];
		dest[2]=base0[2];
		dest[3]=base2[0];
		dest[4]=base2[1];
		dest[5]=base2[2];
		dest[6]=base1[0];
		dest[7]=base1[1];
		dest[8]=base1[2];

	};

	static void mult(VECTOR dest,VECTOR src,double transform[9]){
		for(int i=0;i<3;i++){
			dest[i]=0;
			for(int j=0;j<3;j++)
				dest[i]+=transform[i*3+j]*src[j];
		}
	};
		
	static bool cross(double& result,LINE line0,LINE line1,const double eps){
	//line0: line
	//line1: line segment
		VECTOR dir[2]={
			{
				line0[1][0]-line0[0][0],
				line0[1][1]-line0[0][1],
				line0[1][2]-line0[0][2]
			},{
				line1[1][0]-line1[0][0],
				line1[1][1]-line1[0][1],
				line1[1][2]-line1[0][2]
			}
		};

		VECTOR normal;
		outerProduct(normal,dir[0],dir[1]);

		double transform[9];
		if(!getTransform(transform,dir[0],dir[1],normal))
			return false;

		LINE line0_,line1_;
		for(int i=0;i<2;i++){
			mult(line0_[i],line0[i],transform);
			mult(line1_[i],line1[i],transform);
		}

		if(fabs(line0_[0][2]-line1_[0][2]) > eps)
			return false;
		
		double t=(line1_[0][0]-line0_[1][0])/(line0_[0][0]-line0_[1][0]);
		double s=(line0_[0][1]-line1_[1][1])/(line1_[0][1]-line1_[1][1]);

		if( s<0 || 1<s )
			return false;
		result=t;
		return true;
	};
};

