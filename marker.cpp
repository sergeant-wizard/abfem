/*
ver.03.03
file marker.cpp
*/

#include "marker.h"

Marker::Marker(void):
weight(1.),
isSurface(false),
parent_elem(0),
material(NULL)
{
	F.identity();
	F_prev.identity();
	F_.zero_out();
	F_prev_.zero_out();
}

void Markers_handler::read_gid(ostringstream& filename){
	using std::ifstream;
	using std::cout;
	using std::endl;
	static const int MAX_CHAR=0x100;
	//surface marker
	ifstream fin(filename.str().c_str());
	//issue5
	try{
		if(fin.fail()){
			cout << "marker file not found" << endl;
			throw;
		}
	}
	catch(int){}
	char dump[MAX_CHAR];
	char title[MAX_CHAR];

	fin.getline(dump,MAX_CHAR);
	//title
	fin >> dump >> title;
	for(int i=0;i<3;i++)
		fin.getline(dump,MAX_CHAR);

	//coordinates
	fin >> dump >> num_markers;
	markers=new Marker[num_markers];
	
	for(int i=0;i<2;i++)
		fin.getline(dump,MAX_CHAR);

	for(int i=0;i<num_markers;i++){
		fin >> dump 
			>> markers[i].x[0]
			>> markers[i].x[1]
			>> markers[i].x[2];
		for(int j=0;j<3;j++)
			markers[i].X[j]=markers[i].x[j];
		fin.getline(dump,MAX_CHAR);
	}

	for(int i=0;i<2;i++)
		fin.getline(dump,MAX_CHAR);

	//border
	for(int i=0;i<2;i++)
		fin.getline(dump,MAX_CHAR);

	while(1){
		Connectivity* newCon=new Connectivity;
		fin >> newCon->elem_index >>
			newCon->node_indices[0] >>
			newCon->node_indices[1] >>
			newCon->node_indices[2];
		fin.getline(dump,MAX_CHAR);
		if(fin.fail()){
			delete newCon;
			fin.clear();
			break;
		}

		for(int i=0;i<NODE_IN_ELEM;i++){
			newCon->node_indices[i]--;
			newCon->pmarker[i]=&markers[newCon->node_indices[i]];
		}
		connectivity.push_back(newCon);
	}

	num_faces=connectivity.size();
	try{
		if(!num_faces){
			cout << "no connectivity found" << endl;
			throw;
		}
	}
	catch(int){}

	boundary=new Face[num_faces];
	for(unsigned int i=0;i<num_faces;i++){
		for(int j=0;j<NODE_IN_ELEM;j++){
			markers[connectivity[i]->node_indices[j]].weight=.5;
			markers[connectivity[i]->node_indices[j]].isSurface=true;
		}
	}
}

Markers_handler::~Markers_handler(void){
	delete[] markers;
	delete[] boundary;
	for(unsigned int i=0;i<connectivity.size();i++)
		delete connectivity[i];
}

void Markers_handler::update_element(void){
	for(unsigned int i=0;i<velements.size();i++)
		velements[i]->marker_clear();

	for(int i=0;i<num_markers;i++){
		if(velements[markers[i].parent_elem]->isInside(markers[i].x)){
			velements[markers[i].parent_elem]->marker_add(&markers[i]);
		}else{
			markers[i].parent_elem=-1;
		}
	}

#ifdef __INTEL_COMPILER
	#pragma omp parallel for
#endif
	for(unsigned int j=0;j<velements.size();j++){
		for(int i=0;i<num_markers;i++){
			if(markers[i].parent_elem>=0)
				continue;
			if(velements[j]->isInside(markers[i].x)){
				markers[i].parent_elem=j;
				velements[j]->marker_add(&markers[i]);
			}
		}
	}

#ifdef CHECK_MARKER
	using std::cout;
	using std::endl;
	try{
	for(int i=0;i<num_markers;i++){
		if(markers[i].parent_elem<0){
			cout<<"error: marker outside background mesh"<<endl;
			cout << "marker" << i << endl;
			for(int j=0;j<3;j++)
				cout << markers[i].x[j] << "\t";
			cout<<endl;
			for(int j=0;j<3;j++)
				cout << markers[i].X[j] <<"\t";
			cout<<endl;
			throw;
		}
	}
	}catch(int){return;}
#endif

#ifdef __INTEL_COMPILER
	#pragma omp parallel for
#endif
	for(unsigned int i=0;i<num_faces;i++)
		for(int j=0;j<3;j++)
			for(int k=0;k<3;k++)
				boundary[i].initialize(
					markers[connectivity[i]->node_indices[j]].x[k]
				);

#ifdef __INTEL_COMPILER
	#pragma omp parallel for
#endif
	for(unsigned int j=0;j<velements.size();j++)
		for(unsigned int i=0;i<num_faces;i++)
			if(velements[j]->isInside(boundary[i]))
				velements[j]->face_add(&boundary[i]);
}

