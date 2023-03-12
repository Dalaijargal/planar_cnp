#include "cndp_greedy1_spr.h"

using namespace std;

cndp_greedy1_spr::cndp_greedy1_spr(int k1, int k2, graph_type &component, CNDP_Graph *graph_obj){
	this->total_k = k1;
	this->budget_k = k2;
	this->thecomponent = component;
	this->graph_object = graph_obj;
	this->dim = graph_obj->dim;
}

set <int> cndp_greedy1_spr::greedy1(){
	int st = clock();
	set <int> S = this->find_vertex_cover();
	int tt = clock();
	//cout<<"\nrunning time for VC: "<< (tt-st)/double(CLOCKS_PER_SEC) <<"s\n\n"<< endl;

	set <int> V, V_S;
	// vector <int> cid(dim, -1); // i-th node in cid[i]-th component

	// find V_S
	for (int i = 0; i < dim; i ++)	V.insert(i);
    std::set_difference( V.begin(), V.end(), S.begin(), S.end(), inserter(V_S, V_S.begin()));

    //cout<<"V_S size: "<<V_S.size()<<endl;
    //cout<<"V size: "<<V.size()<<endl;
    //cout<<"S size: "<<S.size()<<endl;

    // define cid[i]
    separator spr_gr;
	vector <int> vc_tmp;
	vc_tmp.reserve(S.size());
	for (set <int>::iterator it = S.begin(); it != S.end(); it ++){
		vc_tmp.push_back(*it);
	}


	//cout<<"vc_tmp size: "<<vc_tmp.size()<<endl;

//	exit(0);

	spr_gr.vertex_list.push_back(vc_tmp);
	//exit(0);
	solution sln_gr = this->graph_object->rmake_solution(spr_gr);
//	exit(0);
	//this->graph_object->print_comp_info(sln_gr);

	//exit(0);

	vector <int> cid = this->graph_object->find_connected_components2(sln_gr.values);
	int num = *std::max_element(cid.begin(), cid.end()) + 1;

	//
	vector <set<int> > comps_set; // sets of component
	comps_set.resize(100*dim);
	for (int i = 0; i < dim; i++){
		if (cid.at(i) != -1){
			comps_set.at(cid.at(i)).insert(i);
		}
	}
	//cout<<"# of components: "<<num<<endl;
	std::vector<int>::size_type sz, sz2;
	sz = comps_set.capacity();
	sz2 = comps_set.size();



	/*
	 * Main part: remove nodes whose removal causes min increase in cost function value
	 * repeat this step until the number of nodes in S meets requirement
	 */
	bigInt_t fn = 0;
	bigInt_t best_min = INT_MAX;
	// add back v in S to G one-by-one
	set <int>::const_iterator sit = S.begin();
	set <int>::const_iterator tit = S.end();
	int counter = 0;
	int counter_newcomp = 0;
	int counter_reall = 0;
	set<int> incidentC;
	// 1st loop
	while(S.size() > this->budget_k){
		// find a node in S whose removal causes best improvement on cost function value, and add it to V_S
		bigInt_t delta_cost;
		bigInt_t fn_min = INT_MAX;
		int v_id;
		set<int> incidentC_best;

		// 2nd loop: iterate all nodes in S
	    for (set <int>::const_iterator v = S.begin(); v != S.end(); v ++){
	    	// find all incident components when node v is added to G

	    	incidentC.clear();
	    	//cout<<*v<<": "<<endl;
	    	for (list <int>::const_iterator w = thecomponent.at(*v).adjList.begin(); w != thecomponent.at(*v).adjList.end(); w ++){
	    		//cout<<*w<<" "<<cid.at(*w)<<endl;
	    		if (cid.at(*w) != -1)
	    			incidentC.insert(cid.at(*w));
	    	}

	    	//exit(0);
	    	// calculate the cost if node v is removed from S
	    	//cout<<*v<<":	";
	    	if (incidentC.size() == 0){
	    		delta_cost = 0;
	    		//cout<<0<<" "<<delta_cost;
	    	}
	    	else {
		    	int cost1 = 0, cost2 = 0, c_total = 1;
		    	for (set<int>::iterator it = incidentC.begin(); it != incidentC.end(); it ++){
		    		int n = comps_set.at(*it).size();
		    		//cout<<"\n	incidents: "<<*it<<"; # of elements: "<<n;
		    		c_total += n;
		    		cost1 += (n*(n-1))/2;
		    	}
		    	cost2 = (c_total*(c_total-1)) / 2;
		    	delta_cost = cost2 - cost1;
		    	//cout<<"\n	cost2: "<<cost2<<"; cost1: "<<cost1<<"; delta: "<<delta_cost;
	    	}
	    	//cout<<"\n";

	    	if (delta_cost < fn_min){
	    		fn_min = delta_cost;
	    		v_id = *v;
	    		incidentC_best = incidentC;
	    		sit = v;
	    		if (fn_min <= best_min){
	    			break;
	    		}
	    	}

	    	//exit(0);
	    } // end second order for-loop

	    // test for speeding
	    //set<int>::iterator it = S.begin();
	    //v_id = *it;

	    sit ++;
	    counter ++;
	    best_min = fn_min;
	    //cout<<v_id<<": "<<fn_min<<endl;
    	if((sit == tit)||(counter > S.size()/10)){
    		sit = S.begin();
    		//cout << "\n\nRunning time for VC:	 " << (tt-st)/double(CLOCKS_PER_SEC) <<"s\n\n"<< endl;
    		//return S;
    		//exit(0);
    	}

	    // remove the node v_id from S
	    S.erase(v_id);

	    //cout<<*(++sit)<<endl;
	    //exit(0);

	    // update the data and variables
	    // 1. update cid[i] by assigning the incidentC[0] that is the first component ID.
	    // case 1: there is no incident component. it means the node v_id is isolated. it will be considered as a new component
	    if (incidentC_best.size() == 0){
	    	// create a new component at the end of the components' set
	    	set<int> new_comp;
	    	new_comp.insert(v_id);
	    	comps_set.push_back(new_comp);
	    	cid.at(v_id) = comps_set.size() - 1;
	    	counter_newcomp ++;
	    }
	    else{
		    std::set<int>::iterator it1 = incidentC_best.begin();
		    int new_id = *it1;
		    comps_set.at(new_id).insert(v_id);
		    cid[v_id] = new_id;
		    it1 ++;
		    for (std::set<int>::iterator it = it1; it != incidentC_best.end(); it ++){
		    	for (std::set<int>::iterator it2 = comps_set.at(*it).begin(); it2 != comps_set.at(*it).end(); it2 ++){
		    		// update cid[] by new_id
		    		cid[*it2] = new_id;
		    		// merge the components to the first component
		    		comps_set.at(new_id).insert(*it2);
		    	}
		    	// remove all incident components to v_id except the first component in order
		    	comps_set.at(*it).clear();
		    }
	    }
	    if (sz != comps_set.capacity()){
	    	sz = comps_set.capacity();
	    	counter_reall ++;
	    	//cout<<"re-allocation!!!"<<endl;
	    	//exit(0);
	    }
	}

	//cout<<"# of creating new comp: "<<counter_newcomp<<endl;
	//cout<<"# of reall: "<<counter_reall<<endl;

	//cout<<"size: "<<sz2<<"; capacity: "<<sz<<endl;
	//cout<<"size: "<<comps_set.size()<<"; capacity: "<<comps_set.capacity()<<endl;
	//exit(0);
	return S;
}

set <int> cndp_greedy1_spr::find_vertex_cover(){
	//srand (time(0));
	set <int> vc;
	// convert graph_type to Graph_set_t graph
	Graph_set_t boost_g;
	for (int i = 0; i < dim; i ++){
		if (this->thecomponent.at(i).adjList.size() != 0){
			boost::add_vertex(i, boost_g);
		}
	}
    for (int i = 0; i < this->thecomponent.size(); i++){
    	for (list<int>::iterator it = this->thecomponent.at(i).adjList.begin(); it != this->thecomponent.at(i).adjList.end(); it ++){
    		if (i < *it){
   				boost::add_edge(i, *it, boost_g);
    		}
    	}
    }
    //cout<<"boost graph is created";

    // discover approximate vertex cover
	boost::graph_traits<Graph_set_t>::out_edge_iterator ei, ei_end;
	Edge e;
	bool b;
	vc = std::set<int>();
	boost::graph_traits<Graph_set_t>::edge_iterator e1, e1_end;
	int v, u;
	int n;
	while (boost::num_edges(boost_g) > 0){
			// 1. select an edge at random
			boost::tie(e1, e1_end) = boost::edges(boost_g);
			n = boost::num_edges(boost_g);
			int r = (rand()*rand()) % (n);
			//cout<<"\nn: "<<n<<"; r: "<<r;
			while (r > 0){
				e1 ++;
				r --;
			}

			//cout<<":=> "<<*e1;
			v = boost::source(*e1, boost_g);
			u = boost::target(*e1, boost_g);

			// 2. restore the points of the selected edge
			vc.insert(u);
			vc.insert(v);

			// 3. remove all incident edges on either u or v
			boost::clear_vertex(u, boost_g);
			boost::clear_vertex(v, boost_g);
	}
	//cout<<"\n# of nodes in vc: "<<vc.size()<<endl;
	//for (std::set<int>::iterator set_it = vc.begin(); set_it != vc.end(); set_it ++) cout<<*set_it<<"\n";
	return vc;
}


bigInt_t F(int n){
	return (n*(n-1))/2;
}


vector <BFStree> cndp_greedy1_spr::create_bfs_tree_vector(solution sln){

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
