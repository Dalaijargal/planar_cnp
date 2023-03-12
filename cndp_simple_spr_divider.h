#include <list>
#include "cndp_cycle_spr.h"

class cndp_simple_spr_divider{
	int dim;
	int total_k;
	int budget_k;

	graph_type thecomponent;
	CNDP_Graph *graph_object;
	gamma_bound_t beta_bound;

public:
	cndp_simple_spr_divider(int bdgt_k, int avail_k, graph_type  &component, CNDP_Graph *graph_obj, gamma_bound_t bound);//, CNDP_Graph graph_obj);

public:
	separator execute(int Nc, separator &pre_spr); // return cost
	//separator create_nc_divider(graph_type &g, int root, int nc);
	std::vector <BFStree> create_bfs_tree_vector_copy(solution sln);
	//separator preprocessing();
	//separator sqrt_n_divider();
	//std::set <int> find_vertex_cover();
	//std::set <int> greedy1();
	//std::pair<separator, int> execute_balanced(int nc, int gamma); // return cost
	//separator create_nc_divider_balanced(graph_type &g, int root, int nc, int gamma);
};
