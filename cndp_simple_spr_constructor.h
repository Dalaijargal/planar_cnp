#include "cndp_simple_spr_cutter.h"

class cndp_simple_spr_constructor{

	int dim;
	int total_k;
	int budget_k;
	double roo;
	separator spr_priori;
	graph_type thecomponent;
	CNDP_Graph *graph_object;

public:
	cndp_simple_spr_constructor(int total_k, int buddet_k, graph_type  &component, CNDP_Graph *graph_obj, separator &spr_priori, double ro);//, CNDP_Graph graph_obj);

public:
	std::pair<separator, int> execute_nc(graph_type &component);
	//int_range_t find_nc_range(graph_type &component);
	//int_range_t find_nc_range(graph_type &comp);
	//separator execute_balanced(graph_type &component, int nc, int gamma);
	void set_thecomponent(graph_type &component);
};
