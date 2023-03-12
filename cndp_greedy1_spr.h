#include "CNDP_Graph.h"

class cndp_greedy1_spr{

	int dim;
	int total_k;
	int budget_k;
	graph_type thecomponent;
	CNDP_Graph *graph_object;

public:
	cndp_greedy1_spr(int total_k, int buddet_k, graph_type  &component, CNDP_Graph *graph_obj);
	std::set <int> greedy1();
	std::set <int> find_vertex_cover();
	std::vector <BFStree> create_bfs_tree_vector(solution sln);
public:

};

