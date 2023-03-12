#include <algorithm>
#include <string.h>
#include <cmath>
#include <vector>

#include "cndp_simple_spr_divider.h"

using namespace std;
// constructor class
cndp_simple_spr_divider::cndp_simple_spr_divider(int totalk, int remainingk, graph_type &component2divide, CNDP_Graph *graph_obj, gamma_bound_t bound){
	this->dim = graph_obj->dim;
	this->total_k = totalk;
	this->budget_k = remainingk;
	this->thecomponent = component2divide;
	this->graph_object = graph_obj;
	this->beta_bound = bound;
}

// complexity O(V)
separator cndp_simple_spr_divider::execute(int Nc, separator &pre_spr){

	separator answer = pre_spr;

    /*
     * Second part: divide the graph with FCS result
     */
	// find all connected components and create bfs tree for all components
    solution sln_priori = this->graph_object->rmake_solution(answer);
	vector<BFStree> bfst_vector = this->create_bfs_tree_vector_copy(sln_priori);

	//exit(0);
	for (int i = 0; i < bfst_vector.size(); i ++){
		bool isFound  = false;
		int levels = 0;
		int last_level = 0;
		int level_id;
		double weighted_nc = bfst_vector.at(i).weight * (Nc+1);

		// cout<<"weight: "<<bfst_vector.at(i).weight<<"; weighted_nc: "<<weighted_nc<<endl;
		// divide if weighted nc for this component is more than one. otherwise, the component is short enough
		if (weighted_nc > 1){
			int tao = bfst_vector.at(i).height / weighted_nc;

			//cout<<"height: "<<bfst_vector.at(i).height<<"; "<<bfst_vector.at(i).weight<<endl;
			//cout<<"i: "<<i<<"; tao: "<<tao<<"; weighted_nc: "<<weighted_nc<<endl;
			//cout<<"levels: ";
			//exit(0);
			do{
				levels = levels + tao;
				if (levels < bfst_vector.at(i).height-tao/2-2){
					isFound = true;
					level_id = levels;
				}

				if (isFound){
					answer.vertex_list.push_back(bfst_vector.at(i).tree.at(level_id));
					isFound = false;
				}
			}while(levels <bfst_vector.at(i).height-tao/2-2);
		}
	}// end of for
	return answer;
}

vector <BFStree> cndp_simple_spr_divider::create_bfs_tree_vector_copy(solution sln){

	vector<BFStree> tree_vec;
	// extracting connected components based on divider
	//solution sln;
	//sln.values.resize(dim,false);
	vector <int> comp_index = this->graph_object->find_connected_components(sln.values);//(comp, sln_dvr.values);
	vector <set<int> > comp_sets; // sets of component
	comp_sets.resize(dim);
	for (int i = 0; i < dim; i++){
		if (comp_index.at(i) != -1){
			comp_sets.at(comp_index.at(i)).insert(i);
		}
	}



	// create and store bfs tree for components whose size is bigger than 5.
	int num = 0;
	int acc_height = 0;
	for (int i = 0; i < dim; i++){
		if (comp_sets.at(i).size() != 0){
			num ++;
			BFStree tmp_tree;
    		graph_type c;
    		c.resize(dim);
    		// create a graph_type c for i-th component
    		// iterate each nodes in i-th component
    		int num_el = 0;
    		for (set<int>::iterator v = comp_sets.at(i).begin(); v != comp_sets.at(i).end(); v ++){
    			//cout<<*v<<" ";
    			num_el ++;
    			for (list<int>::iterator v1 = thecomponent.at(*v).adjList.begin(); v1 != thecomponent.at(*v).adjList.end(); v1 ++){
    				if (comp_index.at(*v1) == i){
    					c.at(*v).adjList.push_back(*v1);
    				}
    			}
    		}

    		// restore the discovered component to comp_set
    		if (num_el > 100){
    			//cout<<num_el<<" ";
       			tmp_tree.root = this->graph_object->rfind_tallest_tree_root(c);
       			tmp_tree.num_nodes = num_el;
       			tmp_tree.vertex_set = comp_sets.at(i);
       			tmp_tree.tree = this->graph_object->rcreate_bfs_tree(c, tmp_tree.root);
       			tmp_tree.height = tmp_tree.tree.size();
       			tmp_tree.widest = 0;
       			for (int j = 0; j < tmp_tree.tree.size(); j ++){
       				if (tmp_tree.tree.at(j).size() > tmp_tree.widest)
       					tmp_tree.widest = tmp_tree.tree.at(j).size();
       			}
       			tree_vec.push_back(tmp_tree);
       			acc_height += tmp_tree.height;
    		}
		}
	}
	//cout<<endl;
	//calculate weight for tree height
	for (int i = 0; i < tree_vec.size(); i ++){
		//cout<<"h: "<<tree_vec.at(i).height<<"; acc_h: "<<acc_height<<endl;
		tree_vec.at(i).weight = (double)tree_vec.at(i).height / (double)acc_height;
		//cout<<"h: "<<tree_vec.at(i).height<<"; acc_h: "<<acc_height<<"; h/acc_h: "<<tree_vec.at(i).weight<<endl;
		if (tree_vec.at(i).weight == 0)
			tree_vec.at(i).weight = 0.01;
	}
	//tree_vec.at(0).weight = 1;
	//
	//cout<<"# of 1 in sln_fcs: "<<sln.nb_nodes<<endl;
	//cout<<"# of components with fcs: "<<num<<endl;
	//cout<<"tree vec size: "<<tree_vec.size()<<endl;
	//exit(0);
	return tree_vec;
};

/*
// complexity O(V)
separator cndp_simple_spr_divider::create_nc_divider(graph_type &g, int root, int nc){

	separator answer;
	vector < vector<int> > bfs_tree;

	bfs_tree = this->graph_object->rcreate_bfs_tree(g, root);
	int tree_height = bfs_tree.size();
	int step  = tree_height / nc;
	if (step <= 0){
		//cout<<"Error!!! step is not positive."<<endl;
		return answer;
		exit(0);
	}

	//cout<<"h: "<<tree_height<<endl;
	//cout<<"step: "<<step<<endl;
	//* create a beta divider with beta-width

	bool isFound  = false;
	int levels = 0;
	int last_level = 0;
	int level_id;
	// O(sqrt(V))
	do{
		// find a level meeting requirement
		levels = levels + step;

		//cout<<"~~"<<levels<<endl;
		if (levels < tree_height - step/3 - 2){
			isFound = true;
			level_id = levels;
		}

		// restore level_id-th bfs_tree level into divider
		if (isFound){
			last_level 	= level_id;
			isFound 	= false;
			//cout<<last_level<<": "<<bfs_tree.at(last_level).size()<<" := ";
			vector <int> vec;
			vec.reserve(bfs_tree.at(last_level).size());
			//vec.push_back(last_level);
			//vec.push_back(bfs_tree.at(last_level).size());
			for (int t = 0; t < bfs_tree.at(last_level).size(); t ++){
				//cout<<bfs_tree.at(last_level).at(t)<<" ";
				vec.push_back(bfs_tree.at(last_level).at(t));
			}
			answer.vertex_list.push_back(vec);
			//cout<<endl;
		}
	}while((levels < tree_height - 2) );

	return answer;
}
*/
