#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <sstream>

#include "CNDP_Graph.h"

using namespace std;

bigInt_t mytime = 0;
vector<bool> visited;
vector<int> components;

CNDP_Graph::CNDP_Graph(const char *file_name, int nb_critical_nodes){
    this->K 				= nb_critical_nodes;
    this->input_file_name 	= file_name;
    this->totalConnections 	= 0;
    this->dim 				= 0;
}

bigInt_t CNDP_Graph::getTotalConnections(){
	return totalConnections;
}

void CNDP_Graph::initializeTotalConnections(bigInt_t c){
	totalConnections = c;
}

void CNDP_Graph::setTotalConnections(bigInt_t c){
	totalConnections = totalConnections + c;
}

void CNDP_Graph::set_thecomponent(graph_type component){
	this->thecomponent = component;
}
graph_type CNDP_Graph::get_graph_pointer(){
	return graph_adj;
}

graph_type CNDP_Graph::find_largest_components(const vector<bool> &solution){

	//cout<<"\n\n hello from find largest components: \n";
	//cout<<"\n\nlarge component: \n";
	//graph_type comp = this->thecomponent;
	//int counter = 0;
    //for (int i = 0; i < dim; i++){
    //	if (comp.at(i).size() != 0){
    //		cout<<i<<": ";
    //		counter ++;
    //		for (list<int>::iterator v = comp.at(i).begin(); v != comp.at(i).end(); v ++){
    //			cout<<*v<<" ";
    //		}
    //		cout<<endl;
    //	}
    //}
    //cout<<"size of component: "<<counter<<endl;

	graph_type lrg_comp;
	vector <bool> sln = solution;
	vector <int> comp_index = this->find_connected_components(sln);
    int num = *max_element(comp_index.begin(), comp_index.end())+1;

    // compute component sizes
    vector<int> vec_counter(num, 0);
    int id;
    for (int i = 0; i < dim; i++){
    	if (comp_index.at(i) != -1){
    		id = comp_index.at(i);
    		vec_counter.at(id) ++;
    	}
    }
    //exit(0);
    // find largest component's id
    int lrg_size = 0;
    int lrg_id 	 = -1;
    for (int i = 0; i < num; i++)
    	if (lrg_size < vec_counter.at(i)){
    		lrg_size = vec_counter.at(i);
    		lrg_id = i;
    	}
    // find nodes in the largest component
    //cout<<"\n\nlargest component: ";
    // create large component with graph_type
    lrg_comp.resize(dim);
    for (int i = 0; i < dim; i++){
    	if (comp_index.at(i) == lrg_id){
    		//cout<<i<<" ";
    		for (list<int>::iterator v = this->graph_adj.at(i).adjList.begin(); v != this->graph_adj.at(i).adjList.end(); v ++){
    			//cout<<*v<<" ";
    			if (comp_index.at(*v) == lrg_id)
    				lrg_comp.at(i).adjList.push_back(*v);
    		}

    	}
    	//cout<<endl;
    }

    /*
    // print out largest component
    cout<<endl;
    for (int i = 0; i < dim; i++){
    	if (lrg_comp.at(i).size() != 0){
    		cout<<i<<": ";
    		for (list<int>::iterator v = lrg_comp.at(i).begin(); v != lrg_comp.at(i).end(); v ++){
    			cout<<*v<<" ";
    		}
    		cout<<endl;
    	}
    }
    */
    //exit(0);
	return lrg_comp;
}
vector<int> CNDP_Graph::find_connected_components2(vector<bool> &solution){
	mytime = 0;
	vector<int> ans_component;
	ans_component.resize(dim, -1);
//	memcpy(visited, solution, dim*sizeof(bool));
//	for (int i = 0; i < solution.size(); i++) visited[i] = solution[i];
	visited = solution;
	for (int i = 0; i < dim; i++) components.at(i) = -1;

	// mark the nodes that is not in the component
	for (int i = 0; i < dim; i++){
		if (!this->thecomponent.at(i).status){
			visited.at(i) = true;
		}
	}

    for (int v=0; v < dim; v++)
    {
        if (visited.at(v) == false)
        {
            DFS_find_components(v);
            mytime ++;
        }
    }

    ans_component = components;
    //cout<<"mytime: "<<mytime<<endl;
//    memcpy(ans_components, components, dim*sizeof(dim));
//    for (int i = 0; i < dim; i++) ans_component[i] = components[i];
    return ans_component;
}

vector<int> CNDP_Graph::find_connected_components(vector<bool> &solution){
	mytime = 0;
	vector<int> ans_component;
	ans_component.resize(dim, -1);
//	memcpy(visited, solution, dim*sizeof(bool));
//	for (int i = 0; i < solution.size(); i++) visited[i] = solution[i];
	visited = solution;
	for (int i = 0; i < dim; i++) components.at(i) = -1;

	// mark the nodes that is not in the component
	for (int i = 0; i < dim; i++){
		if (this->thecomponent.at(i).adjList.size() == 0){
			visited.at(i) = true;
		}
	}

    for (int v=0; v < dim; v++)
    {
        if (visited.at(v) == false)
        {
            DFS_find_components(v);
            mytime ++;
        }
    }

    ans_component = components;
    //cout<<"mytime: "<<mytime<<endl;
//    memcpy(ans_components, components, dim*sizeof(dim));
//    for (int i = 0; i < dim; i++) ans_component[i] = components[i];
    return ans_component;
}

void CNDP_Graph::DFS_find_components(int v)
{
    visited.at(v) = true;
    components.at(v) = mytime;
    //mytime++;
    list<int>::iterator i;
    for(i = this->thecomponent.at(v).adjList.begin(); i != this->thecomponent.at(v).adjList.end(); ++i){
        if(!visited.at(*i)){
        	DFS_find_components(*i);
        }
    }
}

bigInt_t CNDP_Graph::calculate_cost_function_value(vector<bool> solution)
{
	//exit(0);
	mytime = 0;
	totalConnections = 0;
	visited = solution;
	// mark the nodes that is not in the component
	for (int i = 0; i < dim; i++){
		if (this->thecomponent.at(i).adjList.size() == 0){
			visited.at(i) = true;
		}
	}
	//exit(0);
    for (int v=0; v < dim; v++)
    {
        if (visited.at(v) == false)
        {
            DFS_calc_cost(v);
            totalConnections+=((mytime)*(mytime-1));
            mytime = 0;
        }
    }
    return totalConnections/2;
}

void CNDP_Graph::DFS_calc_cost(int v)
{
//	cout<<v<<endl;
    visited.at(v)= true;
    mytime++;
    list<int>::iterator i;
    for(i = this->thecomponent.at(v).adjList.begin(); i != this->thecomponent.at(v).adjList.end(); ++i){
        if(!visited.at(*i)){
        	DFS_calc_cost(*i);
        }
    }
}

Graph_t CNDP_Graph::make_boost_graph(graph_type &comp){
	Graph_t G;
	for (int i = 0; i < dim; i ++){
		if (comp.at(i).adjList.size() != 0){
			boost::add_vertex(i, G);
		}
	}

    for (int i = 0; i < comp.size(); i++){
    	for (list<int>::iterator it = comp.at(i).adjList.begin(); it != comp.at(i).adjList.end(); it ++){
    		if (i < *it){
   				add_edge(i, *it, G);
    		}
    	}
    }

    return G;
}

void CNDP_Graph::print_graph(){

	list<int>::iterator iter;
	for (int i = 0; i < graph_adj.size(); i++){
		cout<<i<<": ";
		for (iter = graph_adj.at(i).adjList.begin(); iter != graph_adj.at(i).adjList.end(); iter ++){
			cout<<*iter<<" ";
		}
		cout<<endl;
	}
}



void CNDP_Graph::loadFile(){


	ifstream myfile (input_file_name);
	string line;
	string str1 = ":";
	int counterLine = 0;
	int tmpdim;

	if(!myfile.good()){
		cout << "Unable to open file";
		exit(0);
	}

	else{
		getline(myfile, line);
		istringstream is(line);
		is >> tmpdim;
		this->dim = tmpdim;
		graph_adj.reserve(dim);
		while (counterLine<dim) {
			getline(myfile, line);
			// erase the part of :
			size_t p = line.find(str1);
			line.erase(0,p+2);
      	  // Using istringstream to read the line into integers.
			istringstream iss(line);
			int next = 0;
			node_out_edges tmp_list;
			tmp_list.status = true;
			while (iss >> next){
//				addEdge(counterLine, next);
				tmp_list.adjList.push_back(next);
			}
			graph_adj.push_back(tmp_list);// = tmp_list;
			counterLine++;
		}// end of while
//		exit(0);
	}// end else
  myfile.close();

  // initialize global variables with dim
  visited.resize(dim, false);
  components.resize(dim, -1);
}// end of loadFile


// moved here from simple cutter class

// convert spr to solution
solution CNDP_Graph::rmake_solution(separator &spr){
	solution answer;
	answer.values.resize(dim, false);
	for (list< vector <int> >::iterator list_it = spr.vertex_list.begin(); list_it != spr.vertex_list.end(); list_it ++){
		for (vector <int>::iterator vec_it2 = list_it->begin(); vec_it2 != list_it->end(); vec_it2 ++){
			if ((*vec_it2 < 0)||(dim <= *vec_it2)){
				this->print_spr(spr);
				cout<<"Error in make_solution: "<<*vec_it2<<endl;
				exit(0);
			}
			answer.values.at(*vec_it2) = true;
		}
	}

	int counterOne = 0;
	for (int i = 0; i < dim; i++){
		if (answer.values.at(i))
			counterOne++;
	}
	// cout<<"# of 1 in solution: "<<counterOne<<endl;
	answer.nb_nodes = counterOne;
	answer.cost = this->calculate_cost_function_value(answer.values);
	return answer;
}

vector<int> CNDP_Graph::find_farthest_nodeall(graph_type &g, int v){
	vector<int> visited(dim);
	vector<int> source(dim);
	vector<int> dest(dim);
	vector<int> answer_vec;
	int src_size = 0;
	int dst_size = 0;
	int id;
//	visited.clear();
	for (int i = 0; i < dim; i++)
		visited.at(i) = -1;

//	for (int i = 0; i < dim; i++)cout<<i<<": "<<visited[i]<<endl;

	list<int>::iterator iterl;
	vector<int>::iterator iterv1;
	vector<int>::iterator iterv2;
	int my_time = 1;
	int answer  = -1;
	source.at(src_size) = v;
	src_size++;
	visited[v] = my_time;
	int src_size2;
//	cout<<source.size()<<endl;
	do{
//		cout<<"source: ";for (int t = 0; t < src_size; t++)	cout<<source[t]<<" ";cout<<endl;
		dst_size = 0;
		my_time++;
		for (int i = 0; i < src_size; i++){
			id = source.at(i);
			for (iterl = g.at(id).adjList.begin(); iterl != g.at(id).adjList.end(); iterl++){
//				cout<<*iter<<" ";
				if (visited.at(*iterl) == -1){
					dest.at(dst_size) = *iterl;
					visited.at(*iterl) = my_time;
					answer = *iterl;
					dst_size++;
					src_size2 = src_size;
//					cout<<*iter<<" ";
				}
			}
//			cout<<endl;
		}
//		cout<<"dest: ";	for (int t = 0; t < dst_size; t++)	cout<<dest[t]<<" ";	cout<<endl;
		for (int k = 0; k < dst_size; k++)
			source.at(k) = dest.at(k);
		src_size = dst_size;
	}while(src_size > 0);

	for (int i = 0; i < src_size2; i++)
		answer_vec.push_back(source.at(i));
	return answer_vec;
}

// return a node from the set of farthest nodes at random
int CNDP_Graph::rfind_tallest_tree_root(graph_type &g){
	int r;
	int counter = 0;
	for (int i = 0; i < dim; i++){
		if (g.at(i).adjList.size() != 0)
			counter ++;
	}
	//cout<<"size of component for finding root: "<<counter<<endl;
	if (counter > 0){
		do{
			r = (rand()*rand()) % dim;
			//cout<<r<<" :"<<g.at(r).size()<<endl;
		}while (g.at(r).adjList.size() == 0);
		r = *g.at(r).adjList.begin();
		vector <int> root = this->find_farthest_nodeall(g, r);

		r = rand() % root.size();
		root = this->find_farthest_nodeall(g, root.at(r));
		r = rand() % root.size();
		return root.at(r);
	}
	else{
		cout<<"Error in tallest tree root finder: component size is less than 1!"<<endl;
		exit(0);
	}
	//r = 0;
	//cout<<"Selected node: "<<root.at(r)<<endl;

}

// complexity: O(V + E)
vector< vector <int> > CNDP_Graph::rcreate_bfs_tree(graph_type &g, int root){
	vector< vector <int> > ans;
	vector <int> answer;
	answer.resize(dim, -1);
	vector <int> tmpLevel;
	vector <int> rootLevel;
	int counter = 0;
	int level = 1;
	vector<int>::iterator it2;
	list<int>::iterator it;

	answer.at(root) = level;
	rootLevel.push_back(-123);
	ans.push_back(rootLevel);
	rootLevel = vector<int>();

	rootLevel.push_back(root);
	level++;
	counter++;

	while (!rootLevel.empty()){
		ans.push_back(rootLevel);
		tmpLevel = vector<int>();
		//		tmpLevel.clear();
		for (it2 = rootLevel.begin(); it2!=rootLevel.end(); it2++){
			for (it = g.at(*it2).adjList.begin(); it!=g.at(*it2).adjList.end(); it++){
				if (answer[*it] == -1){
					answer[*it] = level;
					tmpLevel.push_back(*it);
				}
			}
		}
		rootLevel = vector<int>();
		//		rootLevel.clear();
		for (it2 = tmpLevel.begin(); it2 != tmpLevel.end(); it2++){
			rootLevel.push_back(*it2);
		}
		level++;
	}

	rootLevel.push_back(-321);
	ans.push_back(rootLevel);
	return ans;
}

tree_property CNDP_Graph::rfind_tree_properties(graph_type &g, int root){
	tree_property tp;
	int bound = (int)floor(sqrt(g.size()));
	vector<vector<int> > bfs_tree = rcreate_bfs_tree(g, root);
	int tree_height = bfs_tree.size()-1;// = * max_element(levels.begin(), levels.end())+1;

	int totalCuts = 0;
	int counterCuts = 0;
	int tmp_tw = 0;
	tp.smallest_tree_width = INT_MAX;
	tp.largest_tree_width = 0;
	// find widest and narrowest level
	for (int i = 1; i < tree_height; i++){
		tmp_tw = (int)bfs_tree[i].size();
		if ((tmp_tw) < ((int)bound)){
			totalCuts = totalCuts + bfs_tree[i].size();
			counterCuts ++;
		}
		if (tmp_tw < tp.smallest_tree_width)
			tp.smallest_tree_width = tmp_tw;
		if (tmp_tw > tp.largest_tree_width)
			tp.largest_tree_width = tmp_tw;
	}

	if (counterCuts == 0)
		tp.avg_tree_width = -1;
	else
		tp.avg_tree_width = (double)((double)totalCuts / (double)counterCuts);
	return tp;
}

graph_type CNDP_Graph::make_graph_type(std::set <int> vertex_set){
	graph_type g;
	g.resize(dim);
	// create a graph_type g for input vertex set
	// iterate each nodes in the set
	//int num_el = 0;
	for (set<int>::iterator v = vertex_set.begin(); v != vertex_set.end(); v ++){
		for (list<int>::iterator v1 = graph_adj.at(*v).adjList.begin(); v1 != graph_adj.at(*v).adjList.end(); v1 ++){
			if (vertex_set.find(*v1) != vertex_set.end()){
				g.at(*v).adjList.push_back(*v1);
		//		num_el ++;
			}
		}
	}
	return g;
}
/*
void CNDP_Graph::rclear_separator(separator &p){

	for (int i = 0; i < p.cut_size; i++){
//		p.nodes[i].clear();
		p.vertex_list.clear();
	}
	p.cut_levels_list.clear();
	p.root = -1;
	p.cut_size = 0;
	p.total_nb_nodes = 0;

}
*/
/*
separator CNDP_Graph::copy_spr(separator &source){
	list<vector<int> >::iterator list_vec_it;
	vector<int>::iterator vec_it;
	separator destination;
	destination.vertex_list.clear();
	destination.vertex_list = source.vertex_list;
	//destination.cut_levels_list.clear();
	//destination.cut_size = source.cut_size;
	destination.nb_nodes = source.nb_nodes;
	return destination;
}
*/

void CNDP_Graph::print_component(graph_type &g){
	cout<<"\n--------component----------\n";
	for (int i = 0; i < dim; i ++){
		if (g.at(i).adjList.size() != 0){
			for (list<int>::iterator list_it = g.at(i).adjList.begin(); list_it != g.at(i).adjList.end(); list_it ++){
				cout<<*list_it<<" ";
			}
			cout<<endl;
		}
	}
	cout<<"----------------------------\n";
}

void CNDP_Graph::print_comp_info(solution &sln){
//		cout<<"---------------------------------"<<endl;
	    vector<int> my_component_index;
	    vector<int>::iterator vec_it;
	    //my_component_index = this->graph_object->find_connected_components(sln.values);// find_connected_components(sln.values);
	    my_component_index = this->find_connected_components(sln.values);//(this->graph_adj, sln.values);

	    //	    for (vec_it = my_component_index.begin(); vec_it != my_component_index.end(); vec_it++){
	//	    	cout<<*vec_it<<endl;
	//	    }
	    int num = *max_element(my_component_index.begin(), my_component_index.end())+1;

	    vector<int> vec;
	    int counter;
	    for (int j = 0; j < num; j++){
	    	counter = 0;
	    	for (int i = 0; i < dim; i++){
	    		if (my_component_index.at(i) == j)
	    			counter++;
	    	}
	    	if (counter != 0){
	    		vec.push_back(counter);
	    		cout<<counter<<" ";
	    	}
	    }
	    cout<<endl;
	    cout<<"# of components: "<<num;
	    int biggest  = *max_element(vec.begin(), vec.end());
	    int smallest = *min_element(vec.begin(), vec.end());
	    cout<<"; biggest:  "<<biggest;
	    cout<<"; smallest: "<<smallest;
	    float average = accumulate( vec.begin(), vec.end(), 0.0)/vec.size();
	    cout<<"; average:  " << average
	    		//	    	<<"; gamma optimal: "<<this->gamma_opt
	    		//			<<"; gamma min: "<<this->gamma_smallest
			<<endl;

	    cout<<"the cost of solution: "<<sln.cost<<" :"<<sln.nb_nodes<<endl;

   //	    this->gamma_opt = average;
//	    return vec;
}

void CNDP_Graph::print_spr(separator &spr){
	list< vector <int> >::iterator list_it, list_it1;
	vector <int>::iterator vec_it, vec_it2;
	for (list_it = spr.vertex_list.begin(); list_it != spr.vertex_list.end(); list_it++){
		vec_it2 = list_it->begin();
		//vec_it2 ++;
		//vec_it2 ++;
		for (vec_it = vec_it2; vec_it != list_it->end(); vec_it++){
			cout<<*vec_it<<" ";
		}
		cout<<endl;
	}
}

// find all components after removing the solution, and their size vector will be returned
vector<int> CNDP_Graph::get_comp_info(solution &sln){
	    vector<int> my_component_index;
	    vector<int>::iterator vec_it;
	    my_component_index = this->find_connected_components(sln.values);// find_connected_components(sln.values);
	    int num = *max_element(my_component_index.begin(), my_component_index.end())+1;
	    vector<int> vec;
	    int counter;

	    for (int j = 0; j < num; j++){
	    	counter = 0;
	    	for (int i = 0; i < dim; i++){
	    		if (my_component_index.at(i) == j)
	    			counter++;
	    	}
	    	vec.push_back(counter);
//	    	cout<<counter<<" ";
	    }
	    return vec;
}
