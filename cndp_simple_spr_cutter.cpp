#include <algorithm>
#include <string.h>
#include <cmath>
#include <vector>

#include "cndp_simple_spr_cutter.h"

using namespace std;

// constructor class
cndp_simple_spr_cutter::cndp_simple_spr_cutter(int totalk, int remainingk, graph_type &component2divide, CNDP_Graph *graph_obj, double ro){
	this->dim 			= graph_obj->dim;
	this->total_k 		= totalk;
	this->budget_k 		= remainingk;
	this->roo 			= ro;
	this->thecomponent 	= component2divide;
	this->graph_object 	= graph_obj;

	this->gamma_opt 	= -1;
	this->gamma_smallest = INT_MAX;
}

// complexity: O(VlogV) + O(V + E)
separator cndp_simple_spr_cutter::execute(separator divider){//


	/*
	 * Preparation part:
	 * Extract connected components using divider;
	 * Create sub-graphs for each component;
	 * restore then on comp_sets
	 */

	solution sln_dvr = this->graph_object->rmake_solution(divider);//this->rmake_solution(divider);

	vector <bfs_tree_t> comp_set_bfst;
	comp_set_bfst.reserve(100*dim);
	//vector <graph_type> comp_set;
	//comp_set.reserve(2*dim);
	separator answer;

	// extracting connected components based on divider
	vector <int> comp_index = this->graph_object->find_connected_components(sln_dvr.values);//(comp, sln_dvr.values);
	vector <set<int> > comp_sets; // sets of component
	comp_sets.resize(dim);
	for (int i = 0; i < dim; i++){
		if (comp_index.at(i) != -1){
			comp_sets.at(comp_index.at(i)).insert(i);
		}
	}

	// create sub-graphs for each component: Complexity O(nc*V + E) + ?
	double total_atw = 0;
	int Nbig_comp = 0;
    for (int i = 0; i < dim; i++){
    	if (comp_sets.at(i).size() != 0){
    		//cout<<"\n---------------------\n";
    		graph_type c;
    		c.resize(dim);
    		// create a graph_type c for i-th component
    		// iterate each nodes in i-th component
    		int num_el = 0;
    		for (set<int>::iterator v = comp_sets.at(i).begin(); v != comp_sets.at(i).end(); v ++){
    			//cout<<*v<<" ";
    			for (list<int>::iterator v1 = thecomponent.at(*v).adjList.begin(); v1 != thecomponent.at(*v).adjList.end(); v1 ++){
    				if (comp_index.at(*v1) == i){
    					c.at(*v).adjList.push_back(*v1);
    					num_el ++;
    				}
    			}
    		}
    		//cout<<endl;
    		// restore the discovered component to comp_set
    		if (num_el != 0){
    			int root  = this->graph_object->rfind_tallest_tree_root(c);
    			bfs_tree_t bfs_tree = this->graph_object->rcreate_bfs_tree(c, root);
    			comp_set_bfst.push_back(bfs_tree);
        		//comp_set.push_back(c);
        		//cout<<"..\n";
        		if (num_el > 5){
        			Nbig_comp++;
        			//int root  = this->graph_object->rfind_tallest_tree_root(c);
        			tree_property tp = this->graph_object->rfind_tree_properties(c, root);
        			total_atw = total_atw + tp.avg_tree_width;
        			//cout<<"avg_tree_width: "<<tp.avg_tree_width<<endl;
        		}
    		}
    	}
    	//cout<<endl;
    }

	// define gamma limit [gamma.first, gamma.second]
	gamma_bound_t gamma_limit;
	//cout<<"\n\nNbig_comp: "<<Nbig_comp<<" ? "<<comp_set.size()<<"\n\n";
	double atw = total_atw / Nbig_comp;
	double atw_big = 3*atw + total_k/100;
	int Ncomp_big = floor(this->budget_k / atw_big) + 1 + Nbig_comp - 1;
	gamma_limit.second = (dim - this->total_k) / Ncomp_big + this->total_k;

	//cout<<"atw_big: "<<atw_big<<endl;
	//cout<<"Ncomp_big: "<<Ncomp_big<<endl;

	double atw_small = atw - 2*atw / 3 - 2;
	if (atw_small < 1)
		atw_small = 1;
	int Ncomp_small = floor(this->budget_k / atw_small) + 1 + Nbig_comp - 1;
	if (Ncomp_small < 1)
		Ncomp_small = 1;

	//cout<<"atw_small: "<<atw_small<<endl;
	//cout<<"Ncomp_small: "<<Ncomp_small<<endl;

	gamma_limit.first = (dim - this->total_k) / Ncomp_small;

	//cout<<"num of components:	 "<<comp_set.size()<<endl;
	//cout<<"atw: 			"<<atw<<endl;
	//cout<<"gamma limit:	 	["<<gamma_limit.first<<", "<<gamma_limit.second<<"]"<<endl;

	/*
	 * Main part:
	 * call new gamma evaluator for component sets
	 * find gamma creating good cutter under k requirement by examining gamma values with binary search
	 */

	// check gamma_limit
	if (gamma_limit.first > gamma_limit.second){
		cerr<<"gamma_limit is wrong!";
		exit(0);
	}
	// O(VlogV)
	separator spr = this->new_gamma_evaluator(comp_set_bfst, divider, gamma_limit, atw);

return spr;
}

// complexity: O(VlogV)
separator cndp_simple_spr_cutter::new_gamma_evaluator(vector<bfs_tree_t> &comp_set_bfst, separator spr_divider, gamma_bound_t gamma_limit, double atw){

	separator spr_star = spr_divider;
	gamma_bound_t gamma_best;
	double a = 1;//gamma_limit.first;
	double b = this->total_k;//gamma_limit.second;
	//cout<<"a: "<<a<<"; b: "<<b<<endl;
	double m;
	double ro = this->roo;
	double ro2 = 0.5;
	//cout<<"roo: "<<this->roo<<endl;
	int num = 0;
	// complexcity: O(logV*nc*V) = O(VlogV)
	do{
		//cout<<"\n-----------------------\n";
		m = (a + b) / 2;
		gamma_bound_t gm;
		gm.second = m;
		gm.first = ro2*gm.second - ro * atw;
		if (gm.first < 1)
			gm.first = 1;
		//cout<<"["<<gm.first<<", "<<gm.second<<"]; "<<ro<<endl;

		separator spr_acc = spr_divider;
		// complexity: O(nc * V)
		for (int i = 0; i < comp_set_bfst.size(); i ++){

				bfs_tree_t bfs_tree = comp_set_bfst.at(i);

				// complexity: O(V)
				separator spr_result = this->create_gamma_BC(bfs_tree, gm);

				// add new spr to accumulator
				for (list<vector <int> >::iterator list_it = spr_result.vertex_list.begin(); list_it != spr_result.vertex_list.end(); list_it ++){
					spr_acc.vertex_list.push_back(*list_it);
				}
		}

		solution sln_new = this->graph_object->rmake_solution(spr_acc);
		//cout<<"divider + cutter cost: "<<sln_new.cost<<"; "<<sln_new.nb_nodes<<":<=:"<<this->total_k<<endl<<endl;
		//this->graph_object->print_comp_info(sln_new);
		//cout<<"total_k: "<<this->total_k<<endl;
		if (sln_new.nb_nodes <= this->total_k){
			solution sln_best = this->graph_object->rmake_solution(spr_star);
			if (sln_best.cost > sln_new.cost){
				spr_star = spr_acc;
				gamma_best.first = a;
				gamma_best.second = b;
			}
			b = m;
		}
		else
			a = m;
		num++;
	}while((b-a) > 0.5);
	//cout<<"\n-----------------------\n";

	//cout<<"num gamma evaluation: 	"<<num<<endl;
	//cout<<"best gamma: 		"<<gamma_best.first<<" - "<<gamma_best.second<<endl;
	//solution sln = this->graph_object->rmake_solution(spr_star);
	//cout<<"cost: 			"<<sln.cost<<": "<<sln.nb_nodes<<endl;
	//exit(0);
return spr_star;
}

// complexity: O(tree_height + V) = O(V)
separator cndp_simple_spr_cutter::create_gamma_BC(bfs_tree_t &bfs_tree, gamma_bound_t gamma){

	separator answer;
	//this->graph_object->rclear_separator(answer);

	/*
	 * preparation part:
	 * create bfs tree;
	 * define acc vector
	 */
	//vector < vector<int>> bfs_tree;
	vector<int> acc;
	acc.resize(dim+1, 0);
	//bfs_tree = this->graph_object->rcreate_bfs_tree(g, root);
	int tree_height = bfs_tree.size();
	list<int>::iterator iter;
	acc.at(1) = 1;
	for (int i = 2; i < tree_height-1; i++){
		acc.at(i) = 0;
		acc.at(i) = acc.at(i-1) + bfs_tree.at(i).size();
	}

	//for (int i = 0; i < bfs_tree.size(); i++){
	//	for (int j = 0; j < bfs_tree.at(i).size(); j++){
	//		cout<<bfs_tree.at(i).at(j)<<" ";
	//	}
	//	cout<<endl;
	//}
	//exit(0);
	//int count = 0;
	//for (int i = 0 ; i < g.size(); i ++){
	//	if (g.at(i).adjList.size() != 0){
	//		count ++;
	//	}
	//}
	//cout<<"---size of component: "<<count<<endl;

	/*
	 * main part:
	 * find cuts that cut sub-component with gamma sized
	 */
	//cout<<".\n";
	int delta_acc = 0;

	//int comp_size;
	int g_max;// = this->rvns_gamma_max - 2*ceil(tp.avg_tree_width);
	int g_min;
	int l;
	int last_deleted_level;
	int fn = 0;
	int l_min = 0;
	int l_min_id;
	bool isMeet = false;

	// * main loop for create BC
	g_max = gamma.second;
	g_min = gamma.first;
	l = 2;
	last_deleted_level = 0;

	//cout<<g_max<<" ~ "<<g_min<<" "<<endl;
	//cout<<"..\n";
	do{
		delta_acc = 0;
		l_min = INT_MAX;
		isMeet = false;
		// find a level that divides graph with balanced using small-sized cut.
		while ((delta_acc < g_max)&&(l < tree_height-1)){

			delta_acc = acc.at(l-1) - acc.at(last_deleted_level);
			int rest_part = acc.at(tree_height-2) - acc.at(l);
			if ( (g_min <= delta_acc) && (delta_acc <= g_max) && (rest_part > g_min) ){
				//cout<<"rest part: "<<<<endl;
				//fn = ;
				if (l_min >= bfs_tree[l].size()){
					isMeet = true;
					l_min_id = l;
					l_min = bfs_tree[l].size();
					//comp_size = delta_acc;
				}
			}
			l ++;
		}

		//restore it to answer as a cutter if the new discovered level meet the requirement,
		if ((isMeet) ){//&& (answer.total_nb_nodes + bfs_tree[cut_id].size() <= this->budget_k )

			// restore the new discovered level
			l = l_min_id+2;
			vector <int> vec;

			vector<int>::iterator it_vec;
			for (it_vec = bfs_tree.at(l_min_id).begin(); it_vec != bfs_tree.at(l_min_id).end(); it_vec++){
				vec.push_back(*it_vec);
			}
			answer.vertex_list.push_back(vec);
			last_deleted_level = l_min_id;
			answer.nb_nodes += bfs_tree[last_deleted_level].size();

		}
		l ++;
	}while((l < tree_height-1));
	//cout<<"...\n";
	// add last level as an abstact cut
	//vector <int> v1;
	//v1.push_back(tree_height-1);
	//comp_size = dim - acc.at(last_deleted_level);
	//v1.push_back(comp_size);
	//answer.vertex_list.push_back(v1); // adding the last level
	//answer.cut_size++;


	return answer;
}
