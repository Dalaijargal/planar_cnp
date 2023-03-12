//============================================================================
// Project Name        : ps_cnp_march28_exp01
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
// Feature		:It can work for not connected graphs
//============================================================================

#include <iostream>
#include <ctime>
#include "cndp_simple_spr_constructor.h"

using namespace std;

//const char *infile_name = "E:\\planar_cndp\\dataset_collection\\dataset6_dt\\dt10000.txt";
const char *infile_name = "E:\\planar_cndp\\dataset_collection\\dataset6_dt\\dt1000.txt";
const char *outfile_name = "result.txt";

int K 	= 100;
int dim = 0;

separator preprocessing(int total_k, int budget_k, graph_type  &component, CNDP_Graph *graph_obj);

//spr_priori = preprocessing(K, K, component2divide, &graph_obj);

int main(int argc, char* argv[]){

	/*
	 *  Check the number of parameters
	 */
    if (argc < 2) {
    	cout<<"There is no input parameter. Program will use default parameters: "<<endl;
    }
    string str_name1;
    string str_name2;

	if (argc >= 2){
		str_name1 = argv[1];
		infile_name = str_name1.c_str();
		K = atoi(argv[2]);
		str_name2 = argv[3];
		outfile_name = str_name2.c_str();
		//time_limit = atoi(argv[3]);
	}
	cout <<"input file: 	 "<<infile_name<<" " << endl;
	cout <<"k-critical node: "<<K<<" " << endl;
	cout <<"output file: 	 "<<outfile_name<<" " << endl;

	//system("pause");
	/*
	 * preparation part:
	 */
	srand (time(0));
	CNDP_Graph graph_obj(infile_name, K);
	graph_obj.loadFile();
	dim = graph_obj.dim;

	vector<node_out_edges> g(dim);
	g = graph_obj.get_graph_pointer();
	graph_type component2divide = g;
	graph_obj.set_thecomponent(component2divide);

	int start_time = clock();
	/*
	 * 1. Pre-processing
	 */
	separator spr_priori;
	spr_priori = preprocessing(K, K, component2divide, &graph_obj);
	int stop_time1  = clock();
	cout << "\n***running time for preprocessing: " << (stop_time1-start_time)/double(CLOCKS_PER_SEC) <<"s\n\n"<< endl;

	/*
	 * 2. Construction part
	 */
	// find the remaining K
	solution sln_tmp = graph_obj.rmake_solution(spr_priori);
	int budget_k = K-sln_tmp.nb_nodes;

	cout<<"preprocess K: "<<sln_tmp.nb_nodes<<"; remaining: "<<budget_k<<endl;
	pair<separator, int> constr_result;

	if (budget_k > 0){
		double ro = 2;
		cndp_simple_spr_constructor simple_constructor(K, budget_k, component2divide, &graph_obj, spr_priori, ro);
		constr_result = simple_constructor.execute_nc(component2divide);
	}
	else{
		constr_result.first = spr_priori;
		constr_result.second = -1;
		cout<<"~~~ there are no remaining K for constructor"<<endl;
	}
	int stop_time2  = clock();


	/*
	 * 3. Printing the result
	 */



	solution sln_star = graph_obj.rmake_solution(constr_result.first);
	cout<<"\n******************************\nFinal component info: \n";
	cout<<"K:= 		"<<sln_star.nb_nodes<<endl;
	cout<<"Cost: 		"<<sln_star.cost<<endl;
	cout<<"best nc: 	"<<constr_result.second<<endl;
	graph_obj.print_comp_info(sln_star);
	//cout<<"\n******************************\nInitial spr info: \n";

	cout << "\n\n***Total running time: " << (stop_time2-start_time)/double(CLOCKS_PER_SEC) <<"s\n\n"<< endl;


	std::ofstream outfile;
	outfile.open(outfile_name, std::ios_base::app);

	outfile<<sln_star.cost<<" "<<sln_star.nb_nodes<<" "<<(stop_time2-start_time)/double(CLOCKS_PER_SEC)<<endl;
	outfile.close();
	//system("pause");
	return 0;
}

separator preprocessing(int total_k, int budget_k, graph_type  &comp, CNDP_Graph *graph_obj){
	//cout<<"============= Pre-processing: 1st part==============\n";
	separator answer_spr;
	/*
	 * first part:
	 * 		examine if it is sparse enough;
	 * 		call Greedy1.
	 */


	Graph_t G = graph_obj->make_boost_graph(comp);
	std::vector<int> art_points;
	articulation_points(G, std::back_inserter(art_points));
	//std::cout << "\n\n***num APs: 			" << art_points.size();
	cndp_greedy1_spr greedy_spr(total_k, budget_k, comp, graph_obj);
	if (art_points.size() > dim/10){
		int st = clock();
		set<int> vc = greedy_spr.greedy1();//greedy1();
		int tr = clock();
		//cout << "\n\nRunning time for greedy1():	 " << (tr-st)/double(CLOCKS_PER_SEC) <<"s\n\n"<< endl;

		// calculate the cost
		vector <int> vc_tmp;
		vc_tmp.reserve(vc.size());
		for (set <int>::iterator it = vc.begin(); it != vc.end(); it ++){
			vc_tmp.push_back(*it);
		}
		answer_spr.vertex_list.push_back(vc_tmp);
		//solution sln_greedy1 = graph_obj->rmake_solution(answer_spr);
		//cout<<"cost: "<<sln_greedy1.cost<<"; "<<sln_greedy1.nb_nodes<<endl;
		//graph_obj->print_comp_info(sln_greedy1);
	}
	else {
		cout<<"\ndo not need to call AP-divider!!!"<<endl;
	}


	/*
	 * second part:
	 * 		check whether FCS is needed to call;
	 * 		call FCS if it is needed.
	 */

	clock_t t = clock();

	// create object of cndp cycle spr class
    graph_type component = comp;
    bound_t bound;
    cndp_cycle_spr cycle_spr(total_k, budget_k, bound, component, graph_obj);

    // create bfs tree vector
	solution sln_tmp2 = graph_obj->rmake_solution(answer_spr);
	vector<BFStree> bfst_vect = greedy_spr.create_bfs_tree_vector(sln_tmp2);

	// main loop
	for (int i = 0; i < bfst_vect.size(); i ++){
		// check the condition for calling FCS
		//cout<<"~~~check condition for FCS: "<<i<<": "<<bfst_vect.at(i).widest<<" >? "<<2*sqrt(2*bfst_vect.at(i).num_nodes)
		//	<<" num_nodes: "<<bfst_vect.at(i).num_nodes<<endl;
		if ( (bfst_vect.at(i).widest > sqrt(2*bfst_vect.at(i).num_nodes))
			&&(bfst_vect.at(i).num_nodes > 10*sqrt(dim) ) ){

			// call FCS
			cout<<"\n\n***Do need to call FCS-divider!!!"<<endl;
			graph_type c = graph_obj->make_graph_type(bfst_vect.at(i).vertex_set);
			cycle_spr.set_thecomponent(c);
	    	vector  < Vertex > cycle_spr_result = cycle_spr.execute();

	    	// convert boost to int type
	        vector <int> csr(cycle_spr_result.size());
	        for (int i = 0; i < cycle_spr_result.size(); i ++){
	        	csr.at(i) = (int)cycle_spr_result.at(i);
	        }
	        if (csr.size() != 0){
	        	answer_spr.vertex_list.push_back(csr);
	        }
		}
		else{
			//cout<<"\n\n***Do not need to call FCS-divider!!!"<<endl;
		}
	}// end-for


    //cout<<"\n\n";

    solution sln_tmp = graph_obj->rmake_solution(answer_spr);
    //cout<<"\n# of 1 in sln_fcs: "<<sln_tmp.nb_nodes<<"; cost: "<<sln_tmp.cost<<endl;
    //cout<<"============= End pre-processing ==============\n\n";
//exit(0);
	return answer_spr;

}

