/*
ver.03.03
file face.cpp
*/

#include "face.h"

const double Face::m_eps=1.0E-5;

Face::Face(VECTOR normal,VECTOR point):input_cnt(9){
	for(int i=0;i<3;i++){
		this->normal[i]=normal[i];
		this->point[i]=point[i];
	}
}

void Face::initialize(double src){
	vertice[input_cnt/3][input_cnt%3]=src;
	input_cnt++;

	if(input_cnt<9)
		return;
	input_cnt=0;

	for(int i=0;i<3;i++)
		point[i]=vertice[0][i];

	getNormal(normal,vertice[0],vertice[1],vertice[2]);

	VECTOR base[2]={
		{
			vertice[1][0]-vertice[0][0],
			vertice[1][1]-vertice[0][1],
			vertice[1][2]-vertice[0][2]
		},{
			vertice[2][0]-vertice[0][0],
			vertice[2][1]-vertice[0][1],
			vertice[2][2]-vertice[0][2]
		}
	};
}

double Face::distance(VECTOR src){
	return (
		normal[0]*(src[0]-point[0])+
		normal[1]*(src[1]-point[1])+
		normal[2]*(src[2]-point[2])
	);
}

bool Face::Cross(Face& face1,Face& face2){
	using namespace ABFEM::MATH;
	VECTOR tmp;

	VECTOR dir;
	outerProduct(dir,face1.normal,face2.normal);
	if(fabs(dir[0])+fabs(dir[1])+fabs(dir[2]) < 3.*m_eps)
		return false;

	double transform[9];
	getNormalTransform(transform,face1.vertice[0],face1.vertice[1],face1.vertice[2]);

	VECTOR normal2_,point2_;
	for(int i=0;i<3;i++)
		tmp[i]=face2.point[i]-face1.point[i];

	mult(normal2_,face2.normal,transform);
	mult(point2_,tmp,transform);

	double coef=
		(normal2_[0]*point2_[0]+normal2_[1]*point2_[1]+normal2_[2]*point2_[2])/
		(SQR(normal2_[0])+SQR(normal2_[1]));


	VECTOR point_={coef*normal2_[0],coef*normal2_[1],0};

	double transform_inv[9];
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			transform_inv[i*3+j]=transform[i+j*3];

	VECTOR point;
	mult(tmp,point_,transform_inv);
	for(int i=0;i<3;i++)
		point[i]=tmp[i]+face1.point[i];

	LINE intersection;
	for(int i=0;i<3;i++){
		intersection[0][i]=point[i];
		intersection[1][i]=point[i]+dir[i];
	}


	double t1[2],t2[2];
	if(!face1.cross(t1,intersection))
		return false;
	if(!face2.cross(t2,intersection))
		return false;


	if(MAX(t2[0],t2[1]) - MIN(t1[0],t1[1]) < m_eps)
		return false;
	if(MAX(t1[0],t1[1]) - MIN(t2[0],t2[1]) < m_eps)
		return false;
	return true;
}

bool Face::cross(double result[2],LINE line){
	const static int line_index[3][2]={
		{0,1},{0,2},{1,2}
	};
	LINE line0,line1,line2;
	LINE* pline[3]={&line0,&line1,&line2};
	for(int i=0;i<3;i++)
		for(int j=0;j<2;j++)
			for(int k=0;k<3;k++)
				(*pline[i])[j][k]=vertice[line_index[i][j]][k];

	int cross_cnt=0;
	double t[3];
	for(int i=0;i<3;i++)
		if(cross(t[cross_cnt],line,*pline[i],m_eps))
			cross_cnt++;

	if(cross_cnt<2)
		return false;
	
	double min_t=t[0],max_t=t[0];

	for(int i=1;i<cross_cnt;i++){
		if(t[i]<min_t)
			min_t=t[i];
		if(max_t<t[i])
			max_t=t[i];
	}

	result[0]=min_t;
	result[1]=max_t;
	return true;
}

