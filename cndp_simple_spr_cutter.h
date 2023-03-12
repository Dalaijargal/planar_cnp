#include "cndp_simple_spr_divider.h"

class cndp_simple_spr_cutter{

	int dim;
	int total_k;
	int budget_k;
	double roo;

	double gamma_opt;
	double gamma_smallest;

	graph_type thecomponent;
	CNDP_Graph *graph_object;

public:
	cndp_simple_spr_cutter(int total_k, int remaining_k, graph_type  &component, CNDP_Graph *graph_obj, double ro);//, CNDP_Graph graph_obj);

public:
	separator execute(separator divider); // return cost
	//separator new_gamma_evaluator(std::vector<graph_type> &comp_set, separator spr_divider, gamma_bound_t gamma, double atw);
	separator new_gamma_evaluator(std::vector<bfs_tree_t> &comp_set_bfst, separator spr_divider, gamma_bound_t gamma_limit, double atw);
	//separator create_gamma_BC(graph_type &g, int root, gamma_bound_t gamma);
	separator create_gamma_BC(bfs_tree_t &g, gamma_bound_t gamma);
};
