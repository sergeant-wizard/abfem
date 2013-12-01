/*
ver.03.03
file beta.hpp
*/

class Inc_marker_element_beta:public Inc_marker_element{
public:
	Inc_marker_element_beta(
		int* node_index,
		Node* raw_nodes,
		Default_material* material
	):
		Tetra1(node_index,raw_nodes,material),
		Inc_marker_element(node_index,raw_nodes,material)
	{};

	virtual ~Inc_marker_element_beta(){};

	void update_loadstep(void){
		Tetra1_marker::update_loadstep();
		m_p[0]=0;
	};

	void update_iteration(void){
		if(effective_flag){
			Tetra1::update_F();
			Tetra1::update_dNdx();
			Tetra1_UP::update_iteration_unique();
			update_markers();
		}
	};

	void calcK(Sparse& K){
		if(effective_flag)
			Tetra1_marker::calcK(K);
	};
	void calcK_UP(Sparse& K){
		if(effective_flag)
			Tetra1_UP::calcK_UP(K);
	};
	void calcT(MATRIX& T){
		if(effective_flag)
			Tetra1_marker::calcT(T);
	};
	void calcT_p(MATRIX& T){
		if(effective_flag)
			Tetra1_UP::calcT_p(T);
	};
	void fill_node_flag(vector<unsigned>& node_flag){
		if(!effective_flag)
			return;
		for(unsigned a=0;a<m_num_nodes;a++)
			for(unsigned i=0;i<m_dim;i++)
				node_flag[nodes[a]->dof[i].fmi]=true;
		node_flag[m_pressure_id_]=true;
	};
};

class Global_beta:virtual public Global_inc_marker{
public:
	Global_beta(
		const unsigned loadSteps,
		const double tolerance_equivalence,
		const unsigned max_iteration,
		const unsigned interval,
		const double dirichlet,
		bool symmetric=false
	):
		Global(loadSteps,tolerance_equivalence,max_iteration,interval,dirichlet),
		Global_UP(loadSteps,tolerance_equivalence,max_iteration,interval,dirichlet),
		Global_pstab(loadSteps,tolerance_equivalence,max_iteration,interval,dirichlet),
		Global_marker(loadSteps,tolerance_equivalence,max_iteration,interval,dirichlet),
		Global_inc_marker(loadSteps,tolerance_equivalence,max_iteration,interval,dirichlet)
	{};

	virtual ~Global_beta(){};

	template<class PrePost_type,class Material_type>
	void initialize(
		string projectname,Material_type* front_material,
		const unsigned DIM,const unsigned NODE_IN_ELEM,
		Material_type* back_material,
		string elemface_file
	){
		prepost=new PrePost_type(
			vdof,num_dof_dirichlet,num_dof_nond,T,projectname,velements,
			front_material,3,4,markers_handler,back_material
		);

		Global_UP::initialize_unique(velements);
		Element::Up_cast(Global_marker::velements,velements);

		after_preprocess();
		Global_pstab::initialize_unique<Elem_face_marker>(velements,elemface_file,velem_faces);

		//Element::Up_cast(Global_pstab::velem_faces,velem_faces);
		allocate_K();

		node_effective_flag.resize(num_dof_dirichlet+num_row);
	};

	void update_loadstep(void){
		Global_inc_marker::update_loadstep();

		for(unsigned i=0;i<velements.size();i++)
			velements[i]->effective_flag=false;
		for(unsigned i=0;i<velem_faces.size();i++){
			Tetra1_marker* elem0=velem_faces[i]->get_elem(0);
			Tetra1_marker* elem1=velem_faces[i]->get_elem(1);
			if(elem1){
				if(elem0->get_num_marker() && elem1->get_num_marker()){
					elem0->effective_flag=true;
					elem1->effective_flag=true;
				}
			}
		}

		for(unsigned i=0;i<node_effective_flag.size();i++)
			node_effective_flag[i]=false;
		for(unsigned i=0;i<velements.size();i++)
			velements[i]->fill_node_flag(node_effective_flag);

	};

	void findK(void){
		Global_inc_marker::findK();
		//fixme: make constant
		for(unsigned i=num_dof_dirichlet;i<num_row;i++)
			if(!node_effective_flag[i]){
				unsigned mi=i-num_dof_dirichlet;
				(*K)(mi,mi)=1.0E-5;
			}
	};

protected:
	vector<Inc_marker_element_beta*> velements;
	vector<Elem_face_marker*> velem_faces;
	vector<unsigned> node_effective_flag;
};

