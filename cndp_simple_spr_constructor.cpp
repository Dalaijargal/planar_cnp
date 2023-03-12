#include <algorithm>
#include <string.h>
#include <cmath>
#include <vector>
#include "cndp_simple_spr_constructor.h"

using namespace std;



// constructor class
cndp_simple_spr_constructor::cndp_simple_spr_constructor(int totalk, int remainingk, graph_type &component2divide, CNDP_Graph *graph_obj, separator &spr_prior, double ro){
	this->dim 			= graph_obj->dim;
	this->total_k 		= totalk;
	this->budget_k 		= remainingk;
	this->spr_priori 	= spr_prior;
	this->roo 			= ro;
	this->thecomponent 	= component2divide;
	this->graph_object 	= graph_obj;
	this->graph_object->set_thecomponent(component2divide);
}

void cndp_simple_spr_constructor::set_thecomponent(graph_type &c){
	this->thecomponent = c;
}
//complexity: O(no_imp_iterations*(O(VlogV) + O(V))); no_imp_it = log(sqrt(V))
pair<separator, int> cndp_simple_spr_constructor::execute_nc(graph_type &component){//

	int_range_t nc_range;
	nc_range.first = 0;
	nc_range.second = 10 + this->total_k/10;
	gamma_bound_t bound;
	// create object of simple divider and execute it on input graph with budged_k = K
	cndp_simple_spr_divider simple_divider(this->budget_k, this->budget_k, this->thecomponent, graph_object, bound);

	solution sln_best;
	sln_best.values.resize(dim, false);
	sln_best.cost = ULLONG_MAX;
	separator spr_best;
	int best_nc;

	int num_no_imp = 0;

	//int start_s = clock();
	//cout<<"\n\nRunning ... \n";
	separator spr_divider;
	int Nc = nc_range.first;
	Nc = 5;
	//complexity: O(no_imp_iterations*(O(VlogV) + O(V)))
	//separator spr_pre;// = simple_divider.preprocessing();
	do{
		cout<<"\n\n============================"<<endl;
		/*
		 * 1. divider part
		 */
		cout<<"Nc: 		"<<Nc<<endl;
		// complexity O(V)
		spr_divider = simple_divider.execute(Nc, this->spr_priori);//myvns.rexecute();
		solution sln_divider =  this->graph_object->rmake_solution(spr_divider);
		cout<<"divider.K: 	"<<sln_divider.nb_nodes<<endl;
		cout<<"divider cost:	"<<sln_divider.cost<<endl;
		this->graph_object->print_comp_info(sln_divider);

		// 2. check termination condition
		if (sln_divider.nb_nodes >= this->budget_k){
			cout<<"!!!divider.k exceeds K"<<endl;
			break;
		}

		/*
		 * 3. cutter part based on divider result
		 */
		cout<<"\n------\n";
		int budgetk = this->budget_k - sln_divider.nb_nodes;
		cndp_simple_spr_cutter cutter(this->budget_k, budgetk, this->thecomponent, this->graph_object, this->roo);
		// O(VlogV + E)
		separator simple_cutter_result = cutter.execute(spr_divider);//myvns.rexecute();
		solution sln_final	= this->graph_object->rmake_solution(simple_cutter_result);// cutter.rmake_solution(simple_cutter_result);

		cout<<"divider+cutter.K:	 	"<<sln_final.nb_nodes<<endl;
		cout<<"divider+cutter.Cost:		" <<sln_final.cost<<endl;
//exit(0);

		/*
		 * 4. move
		 */
		if (sln_best.cost > sln_final.cost){
			sln_best = sln_final;
			spr_best = simple_cutter_result;
			best_nc = Nc;
			num_no_imp = 0;
		}
		num_no_imp ++;
		Nc ++;
		//break;
	} while ((num_no_imp < 10)&&(Nc <= nc_range.second));//(Nc < graph_obj.dim); (0);//


	pair<separator, int> result;
	result.first = spr_best; //divider_result.first;//
	result.second = best_nc;
return result;
}
