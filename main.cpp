/*
ver.03.03 -> change stdout
Ryo Miyajima
file main.cpp
Since 2010
*/

/*
	conf file
	0 Gid project name
	1 Material ID
	2 Young Modulus
	3 Poisson's ratio
	4 Load Steps
	5 Tolerance
	6 Max Iteration
	7 Output Interval
	8 Dirichlet
*/
 
/*
	issues
	1		get rid of common.h
	2		calculating residue for incompressibility 
   			with relative volumetric strain
	3		minimum marker weight
	4	D	file dependencies
	5		error handling
	6		non-intel compiler enviroment
	7	D	c_p used in bonet's book
	8	D	do update_marker in class Marker
	9		terminal output handler
	10		updated or total lagrange
	11		capsule for beta
	12		bodyforce
*/

#include <iostream>
#include <fstream>
#include <stdlib.h>

#include <math.h>

#include "incompressible.hpp"
#include "elastic.hpp"
#include "global_pstab.hpp"
#include "gid.h"
#include "gid_marker.hpp"
#include "inc_marker.hpp"
#include "gid.hpp"
#include "hexa27.h"
#include "bodyforce.hpp"

//deprecated
//#include "beta.h"

/*modes
	100 tetra1		compressible
	110 tetra1		incompressible
	210 P1-iso-P2	incompressible
	101 tetra1		compressible	marker
	111 tetra1		incompressible	marker
	6 	incompressible marker (beta)	deprecated
	7 	MODE 2 + Body Force				deprecated
	8 	MODE 3 + Body Force				deprecated
	300 hexa27		compressible
	310 hexa27		incompressible
	301 hexa27		compressible	marker
*/
#define MODE 100

int main(int argc,char* argv[]){
	using namespace std;

	cout << "\x1b[0:0H\x1b[2J";
	const static int MAX_CHAR=0x100;
	const static int NUM_OPTIONS=10;
	cout << "ABFEM ver.03.03" << endl;
	cout << "Last compiled " << __TIME__ << endl;
	if(argc<2){
		cout << "usage: ABFEM conf_file (elemface_file)" << endl;
	}else{
		cout << "Debug...";
#ifndef _DEBUG
		cout << "NO" << endl;
#else
		cout << "YES" << endl;
#endif
		ifstream fin(argv[1]);
		try{
			if(fin.fail())
				throw "bad conf file";
		}catch(string str){
			cout<<str<<endl;
		}
		string buff[NUM_OPTIONS];
		char dummy[MAX_CHAR];
		for(int i=0;i<NUM_OPTIONS;i++){
			fin >> buff[i];
			fin.getline(dummy,MAX_CHAR);
		}
		string problem(buff[0]);
		int 		material_id		=atoi(buff[1].c_str());
		double 		young			=atof(buff[2].c_str());
		double 		poisson			=atof(buff[3].c_str());
		unsigned 	loadsteps		=atoi(buff[4].c_str());
		double 		tolerance		=atof(buff[5].c_str());
		unsigned 	max_iteration	=atoi(buff[6].c_str());
		unsigned 	output_interval	=atoi(buff[7].c_str());
		double 		dirichlet		=atof(buff[8].c_str());

		cout << "Problem:" << problem << endl;

#if MODE==100
		Bonet_comp_ material(young,poisson);
		Global global(
			loadsteps,tolerance,max_iteration,output_interval
		);
		global.initialize<Gid,Tetra1,Elastic<3,MATRIX33> >(
			problem,
			&material,
			3,4,1
		);
#elif MODE==110
		if(argc != 3){
			cout << "specify elemface conf file" << endl;
			return 0;
		}
		Neo_hookean material(young);
		Global_pstab global(
			loadsteps,tolerance,max_iteration,output_interval
		);
		global.initialize<Gid,Incompressible<3,MATRIX33> >(
			problem,
			&material,
			argv[2]//fixme?
		);
#elif MODE==210
		//cube-uap12
		Neo_hookean material(young);
		Global_UP global(
			loadsteps,tolerance,max_iteration,output_interval,1
		);
		global.initialize<Gid,Tetra8,Incompressible<3,MATRIX33> >(
			problem,
			&material,
			3,10,1
		);
#elif MODE==101
		//cube-uap101
		Bonet_comp_ front_material(young,poisson);
		Bonet_comp_ back_material(young*1E-2,poisson);
		Global_marker global(
			loadsteps,tolerance,max_iteration,output_interval,dirichlet
		);
		global.initialize<Gid_marker,Tetra1_marker,Elastic<3,MATRIX33> >(
			problem,
			&front_material,
			3,4,1,
			&back_material
		);
#elif MODE==111
		//incompressible marker
		if(argc != 3){
			cout << "specify elemface conf file" << endl;
			return 0;
		}
		Neo_hookean front_material(young);
		Neo_hookean back_material(young*1E-2);
		Global_inc_marker global(
			loadsteps,tolerance,max_iteration,output_interval,dirichlet
		);
		global.initialize<Gid_marker,Incompressible<3,MATRIX33> >(
			problem,
			&front_material,
			3,4,
			&back_material,
			argv[2]
		);
#elif MODE==6
		//incompressible marker
		if(argc != 3){
			cout << "specify elemface conf file" << endl;
			return 0;
		}
		Neo_hookean front_material(young);
		Neo_hookean back_material(0);
		Global_beta global(
			loadsteps,tolerance,max_iteration,output_interval,dirichlet
		);
		global.initialize<Gid_marker,Incompressible<3,MATRIX33> >(
			problem,
			&front_material,
			3,4,
			&back_material,
			argv[2]
		);
#elif MODE==7
		if(argc != 3){
			cout << "specify elemface conf file" << endl;
			return 0;
		}
		{
			string dump;
			ifstream fin("bodyforce.conf");
			fin >> dump >> g_density;
			fin >> dump >> g_gravity;
			fin >> dump >>
				g_bodyforce_normal[0] >> 
				g_bodyforce_normal[1] >> 
				g_bodyforce_normal[2];
			double size=1./sqrt(
				g_bodyforce_normal[0]*g_bodyforce_normal[0]+
				g_bodyforce_normal[1]*g_bodyforce_normal[1]+
				g_bodyforce_normal[2]*g_bodyforce_normal[2]
			);
			for(unsigned i=0;i<3;i++)
				g_bodyforce_normal[i]*=size;
			g_loadsteps=loadsteps;
		}
		Neo_hookean material(young);
		Global_bodyforce global(
			loadsteps,tolerance,max_iteration,output_interval
		);
		global.initialize<Gid,Incompressible<3,MATRIX33> >(
			problem,
			&material,
			argv[2]//fixme?
		);
#elif MODE==8
		if(argc != 2){
			cout << "specify conf file" << endl;
			return 0;
		}
		{
			string dump;
			ifstream fin("bodyforce.conf");
			fin >> dump >> g_density;
			fin >> dump >> g_gravity;
			fin >> dump >>
				g_bodyforce_normal[0] >> 
				g_bodyforce_normal[1] >> 
				g_bodyforce_normal[2];
			double size=1./sqrt(
				g_bodyforce_normal[0]*g_bodyforce_normal[0]+
				g_bodyforce_normal[1]*g_bodyforce_normal[1]+
				g_bodyforce_normal[2]*g_bodyforce_normal[2]
			);
			for(unsigned i=0;i<3;i++)
				g_bodyforce_normal[i]*=size;
			g_loadsteps=loadsteps;
		}
		Neo_hookean material(young);
		Global_UP global(
			loadsteps,tolerance,max_iteration,output_interval
		);
		global.initialize<Gid,Tetra8,Incompressible<3,MATRIX33> >(
			problem,
			&material,
			3,10
		);
#elif MODE==300
		Bonet_comp_ material(young,poisson);
		Global global(
			loadsteps,tolerance,max_iteration,output_interval
		);
		global.initialize<Gid,Hexa27,Elastic<3,MATRIX33> >(
			problem,
			&material,
			3,27,27
		);
#elif MODE==310
		Neo_hookean material(young);
		Global_UP global(
			loadsteps,tolerance,max_iteration,output_interval,4
		);
		global.initialize<Gid,Hexa27_UP,Incompressible<3,MATRIX33> >(
			problem,
			&material,
			3,27,27
		);
#elif MODE==301
		Venant_ front_material(young,poisson);
		Venant_ back_material (young*1E-3,poisson);

		Global_marker global(
			loadsteps,tolerance,max_iteration,output_interval,dirichlet
		);
		global.initialize<Gid_marker,Hexa27_marker,Elastic<3,MATRIX33> >(
			problem,
			&front_material,
			3,27,27,
			&back_material
		);
#endif

		global.solve();
	}
	cout << "exiting MAFS" << endl;
}

		/*
		Material* material;
		switch(atoi(buff[1].c_str())){
			case 0:
				material=new Bonet1(
					atof(buff[2].c_str()),
					atof(buff[3].c_str())
				);
				break;
			case 1:
				material=new Bonet2(
					atof(buff[2].c_str()),
					atof(buff[3].c_str())
				);
				break;
			case 2:
				material=new Owen(
					atof(buff[2].c_str()),
					atof(buff[3].c_str())
				);
				break;
			case 3:
				material=new Venant_(
					atof(buff[2].c_str()),
					atof(buff[3].c_str())
				);
				break;
			case 4:
				material=new Simo(
					atof(buff[2].c_str())
				);
				break;
		};
		*/

	/*
	double c[3][3][3][3];
	Matrix33 A;
	Matrix33 sigma;
	A(0,0)=-0.0001*0.5;
	A(1,1)=-0.0001*0.5;
	A(2,2)=0.0001;
	Bonet_ material(1.0);

	material.getC_e(c,&A);
	material.getSigma(&sigma,&A);

	sigma.print(cout);

	cout << "A" << endl;
	A.print(cout);

	Matrix33 B;
	B.copy(A);
	cout << "B" << endl;
	B.print(cout);

	cout << c[0][0][0][0] << endl;
	cout << c[0][0][1][1] << endl;
	cout << c[0][1][0][1] << endl;

	return 0;
	*/

	/*
	Node* nodes[4];
	for(unsigned i=0;i<4;i++)
		nodes[i]=new Node(3);
	nodes[0]->dof[0].X=0;
	nodes[0]->dof[1].X=0;
	nodes[0]->dof[2].X=0;

	nodes[1]->dof[0].X=1;
	nodes[1]->dof[1].X=0;
	nodes[1]->dof[2].X=0;

	nodes[2]->dof[0].X=0;
	nodes[2]->dof[1].X=1;
	nodes[2]->dof[2].X=0;

	nodes[3]->dof[0].X=0;
	nodes[3]->dof[1].X=0;
	nodes[3]->dof[2].X=1;

	int indices[4]={1,2,3,4};
	//Venant_ material(1.0,0.3);
	Simo material(1.0);

	//Tetra1_UP elem(indices,*nodes,&material);
	//Tetra1 elem(indices,*nodes,&material);
	//Element* super=&elem;
	//MATRIX_WIZ T(12,1);
	//super->calcT(T);

	Tetra1_UP* elem=new Tetra1_UP(indices,*nodes,&material);
	Element_output* elem_output=elem;

	delete elem_output;
	//delete elem;

	for(unsigned i=0;i<4;i++)
		delete nodes[i];

	return 0;

	*/

/*Mooney-Rivlin test
	using namespace std;
	const double a=1E-4;
	Matrix33 F_;

	F_(0,0)=-a*.5;
	F_(1,1)=-a*.5;
	F_(2,2)=a;

	Matrix33 sigma;
	Mooney_Rivlin mat(1.,1./6.);
	mat.getSigma(&sigma,&F_);

	for(unsigned i=0;i<3;i++){
		for(unsigned j=0;j<3;j++)
			cout << sigma(i,j) << "\t";
		cout << endl;
	}

	double c[3][3][3][3];
	mat.getC_e(c,&F_);
	cout << c[0][0][0][0] << endl;
	cout << c[0][0][1][1] << endl;
	cout << c[0][1][0][1] << endl;

	return 0;
*/
