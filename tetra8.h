/*
ver.03.03
file tetra8.h
*/

#pragma once
#include "tetra1_up.hpp"

class Tetra8:public Element_UP{
public:
	typedef Incompressible<3,MATRIX33> Default_material;
	Tetra8(
		int* indices,
		Node* nodes,
		Default_material* material
	);
	~Tetra8();

	//from Element
	void initialize_sparse(Sparse& K);
	
	void initialize(void){
		for(unsigned i=0;i<ELEM_IN_ELEM;i++)
			subelements[i]->initialize();
	};
	void update_loadstep(void){
		for(unsigned elem=0;elem<ELEM_IN_ELEM;elem++)
			subelements[elem]->update_loadstep();
	};
	void update_iteration(void){
		for(unsigned elem=0;elem<ELEM_IN_ELEM;elem++)
			subelements[elem]->update_iteration();
	};
	void calcK(Sparse& K){
		for(unsigned elem=0;elem<ELEM_IN_ELEM;elem++)
			subelements[elem]->calcK(K);
	};
	void calcT(MATRIX& T){
		for(unsigned elem=0;elem<ELEM_IN_ELEM;elem++)
			subelements[elem]->calcT(T);
	};
	Node*& get_node(unsigned node){
		return nodes[node];
	};

	void sync_Kmatrix(Sparse& K,bool symmetric=true){
		for(unsigned elem=0;elem<ELEM_IN_ELEM;elem++)
			subelements[elem]->sync_Kmatrix(K);
	};

	//from Element_UP
	void put_pressureID(unsigned pi)const{
		Element_UP::put_pressureID(pi);
		for(unsigned elem=0;elem<ELEM_IN_ELEM;elem++)
			subelements[elem]->put_pressureID(pi);
	};
	void put_pressureID_(unsigned pi)const{
		Element_UP::put_pressureID_(pi);
		for(unsigned elem=0;elem<ELEM_IN_ELEM;elem++)
			subelements[elem]->put_pressureID_(pi);
	};
	void calcK_UP(Sparse& K){
		for(unsigned elem=0;elem<ELEM_IN_ELEM;elem++)
			subelements[elem]->calcK_UP(K);
	};
	void calcT_p(MATRIX& T){
		for(unsigned elem=0;elem<ELEM_IN_ELEM;elem++)
			subelements[elem]->calcT_p(T);
	};
	void increment_pressure(unsigned alpha,double dp){
		Element_UP::increment_pressure(alpha,dp);
		for(unsigned elem=0;elem<ELEM_IN_ELEM;elem++)
			subelements[elem]->increment_pressure(alpha,dp);
	};

	//others
	Tetra1_UP* sub(unsigned elem){
		return subelements[elem];
	};

	//issue12
	void calcBodyForce(MATRIX& T);

	static const unsigned ELEM_IN_ELEM=8;

	/*obsolete
	double get_J(void);
	void put_pressureID(unsigned id)const;
	double get_weight(void){
		return volume_ratio;
	};
	void increment_pressure(double dp);
	void calcVolume(void);
	*/
private:
	static const unsigned EDGE[6][2];
	static const unsigned NODE_IN_ELEM=10;
	static const unsigned NODE_SMALL_ELEM=4;
	static const unsigned DIM=3;
	static const unsigned ELEM_CONNECT[ELEM_IN_ELEM][NODE_SMALL_ELEM];
	static const short SUP_TO_SUB[ELEM_IN_ELEM][NODE_IN_ELEM];
	static const double m_eps;
	int pressureID;

	Tetra1_UP* subelements[ELEM_IN_ELEM];
	Default_material* material;
	Node* nodes[NODE_IN_ELEM];

	//obsolete
	double volume_ratio;
};

