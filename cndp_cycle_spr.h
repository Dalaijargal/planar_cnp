#ifndef CNDP_CYCLE_SPR_H
#define CNDP_CYCLE_SPR_H
#include "cndp_greedy1_spr.h"

#include <boost/config.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/boyer_myrvold_planar_test.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/subgraph.hpp>
//#include <boost/graph/planar_detail/planar_dual.hpp>
#include <boost/graph/planar_detail/planar_dual_tree2.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/tuple/tuple.hpp>

#include <boost/graph/make_biconnected_planar.hpp>
#include <boost/graph/make_maximal_planar.hpp>
#include <boost/graph/planar_face_traversal.hpp>
#include <boost/graph/boyer_myrvold_planar_test.hpp>


typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
		boost::property<boost::vertex_discover_time_t , int, boost::property <boost::vertex_index_t, int> >,
		boost::property<boost::edge_weight_t, int, boost::property <boost::edge_index_t, int> > > Graph_t;
typedef typename boost::graph_traits<Graph_t>::edge_descriptor Edge;
typedef typename boost::graph_traits<Graph_t>::vertex_descriptor Vertex;
typedef typename std::set<Edge> EdgeSet;
typedef std::pair<int,int> E;
typedef boost::property_map <Graph_t, boost::edge_weight_t>::type pm_edgePropertyType;
typedef boost::property_traits<pm_edgePropertyType>::value_type pt_edgePropertyType;
typedef boost::graph_traits<Graph_t>::edge_iterator gt_edgeIterator;
typedef boost::property_map<Graph_t, boost::edge_index_t>::type pm_edgeIndexType;
typedef typename boost::property_map<Graph_t, boost::vertex_discover_time_t>::type pm_vertex_disc_t;
typedef std::pair<double, double> bound_t;

struct big_face_t{
	std::vector<bool> table;
	std::vector<int> dist;
	Vertex root;
	int length;
};

struct fsc_info{
	Edge entrance_edge;
	std::list<Vertex> cycle;
	int in;
	int out;
	int spr;
};

struct edge_t{
	std::pair<Edge, int> edge;
};

struct separator_t{
	int out;
	int in;
	int spr;
};

struct dual_tree_node{
	Vertex vertex_id;				// store vertex id in dual graph
	std::vector<edge_t> face_edge_list;	// the corresponding edges of the face
	std::vector<Vertex> face_vertex_list;// store all vertices on the face
	int bfs_level;					// store bfs_level
	Vertex parent_id;				// store parent vertex id on dual tree
	std::vector <Vertex> child_id;		// store child vertex id on dual tree
	Edge entrance_edge;				// store entrance edge in the primal
	std::vector <Edge> outgoing_edges;	// store  outgoing edges in the primal
	separator_t spr_info; 			// separator information for the entrance_edge
	std::set <Vertex> f_cycle;			// store vertices on the fundamental cycle
	Graph_t cycle;
	int cycle_type;
};

class cndp_cycle_spr{
	int total_k;
	int budget_k;
	int dim;
	bound_t bound;
	graph_type thecomponent;
	CNDP_Graph *graph_object;
	//	Graph_t g;
	//	Graph_t g3;
	//	Vertex root;
	//	set<Edge> tree_edges;
	//	set<Edge> tree_edges3;
	//	graph_type component;
	//	CNDP_Graph * graph_obj;
	//	set<Vertex> big_face;

public:
	cndp_cycle_spr(int k1, int k2, bound_t bound, graph_type component, CNDP_Graph *graph_obj);

	std::vector<Vertex> execute();
	std::vector <dual_tree_node> create_dtree();
	std::vector <dual_tree_node>  rank_ntree_edges(Vertex root, Graph_t &g, EdgeSet &my_tree);
	std::vector<dual_tree_node> compute_fcs(Graph_t &g, EdgeSet &tree_edges, Vertex root);

	std::vector <std::vector <Vertex> > find_path_spr(Graph_t &g, std::pair<double, double> bound, big_face_t big_face, EdgeSet &tree_edges);


	Vertex find_o(std::vector<std::pair<Vertex, bool> > &te_tree, std::set<Vertex> &fc_intersect, Vertex p2, Vertex &root);

	Graph_t triangulate(Graph_t &g3);
	std::set<Edge> redifine_tree_edges(std::set<Edge> &tree_edges, Graph_t &g3);
	void find_good_cycles(std::vector <dual_tree_node> &dtree);

	big_face_t find_biggest_face(Graph_t &g);
	int find_height(Graph_t &g, Vertex v);
	Vertex height_minimize(Graph_t &g, std::set<Vertex> &big_face);
	std::vector <std::vector <Vertex> > select_good_cycles(Graph_t &g, std::vector <dual_tree_node> &mytree, std::pair<double, double> bound, big_face_t big_face, EdgeSet &tree_edges);
	std::vector < std::set <Vertex> > clean_fcs_cycle(std::vector < std::set <Vertex> > dtree, std::set<Vertex> &bface, Graph_t &g);
	void set_thecomponent(graph_type &component);

	void print_dual_tree(std::vector<dual_tree_node> tree);
	void print_spr_list(std::list <fsc_info> & my_list);
	template <typename Graph_t>
	void print_graph(Graph_t& g);
	void print_spr_vector(std::vector <std::set <Vertex> > & my_vector);

};

#endif
