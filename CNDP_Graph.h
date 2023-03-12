#include<iostream>
#include <list>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <unordered_set>
#include <string.h>
#include <time.h>
#include <math.h>
#include <algorithm>
#include <numeric>

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
#include <boost/graph/connected_components.hpp>
#include <boost/graph/biconnected_components.hpp>

typedef boost::property<boost::vertex_discover_time_t , int, boost::property <boost::vertex_index_t, int> > VertexProperty;
typedef boost::property<boost::edge_weight_t, int, boost::property <boost::edge_index_t, int> > EdgeProperty;

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, VertexProperty, EdgeProperty > Graph_t;
typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS, VertexProperty, EdgeProperty > Graph_set_t;
typedef typename boost::graph_traits<Graph_t>::edge_descriptor Edge;
typedef typename boost::graph_traits<Graph_t>::vertex_descriptor Vertex;
typedef typename std::set<Edge> EdgeSet;
typedef std::pair<int,int> E;

typedef typename boost::property_map <Graph_t, boost::edge_weight_t>::type pm_edgePropertyType;
typedef typename boost::property_map<Graph_t, boost::edge_index_t>::type pm_edgeIndexType;
typedef typename boost::property_map<Graph_t, boost::vertex_discover_time_t>::type pm_vertexDiscoveredTimeType;

typedef typename boost::property_traits<pm_edgePropertyType>::value_type pt_edgePropertyType;

typedef typename boost::graph_traits<Graph_t>::edge_iterator gt_edgeIterator;
typedef typename boost::graph_traits<Graph_t>::out_edge_iterator gt_outEdgeIterator;

typedef std::pair<double, double> bound_t;

typedef std::pair<double, double> gamma_bound_t;

typedef unsigned long long int bigInt_t;

typedef std::pair<int, int> int_range_t;

typedef std::vector < std::vector<int>> bfs_tree_t;

struct node_out_edges{
	std::list< int > adjList;
	bool status;
};
typedef std::vector< node_out_edges > graph_type;

struct separator{
	int nb_nodes;
	std::list< std::vector <int> > vertex_list;
	std::vector<bool> sln;
};

struct solution{
	bigInt_t cost;
	int nb_nodes;
	std::vector<bool> values;
};

struct component{
	int root;
	std::vector <int> vertex_vector;
	std::list<std::vector <int > > spr_list;
};

struct tree_property{
	double avg_tree_width;
	int smallest_tree_width;
	int largest_tree_width;
};

struct BFStree{
	std::set <int> vertex_set;
	std::vector <std::vector<int> > tree;
	int num_nodes;
	int height;
	int root;
	double weight;
	int widest;
};
// Graph class represents a undirected graph
// using adjacency list representation
class CNDP_Graph
{
    //int V;    // No. of vertices

    int K;
    const char *input_file_name;
    unsigned long long int totalConnections;
    graph_type thecomponent;
    graph_type graph_adj;
   // Pointer to an array containing adjacency lists

   // A function used by DFS
    void DFS_calc_cost(int v); //, bool visited[]
    void DFS_find_components(int v);

public:
    int dim;

public:
    CNDP_Graph(const char *input_file_name, int nb_critical_nodes);   // Constructor
    void addEdge(int v, int w);
    bigInt_t calculate_cost_function_value(std::vector<bool> s);
    std::vector<int> find_connected_components(std::vector<bool> &solution);
    std::vector<int> find_connected_components2(std::vector<bool> &solution);
    graph_type find_largest_components(const std::vector<bool> &sln);
    void loadFile();

    //void deleteData();
//    bool * getBestSolution();
    bigInt_t getTotalConnections();
    graph_type get_graph_pointer();
    void setTotalConnections(bigInt_t c);
    void initializeTotalConnections(bigInt_t c);
    void print_graph();
    void set_thecomponent(graph_type component);

    // move from cutter and divider class
    std::vector<int> get_comp_info(solution &sln);
    void print_spr(separator &spr);
    void print_comp_info(solution &sln);
    void print_component(graph_type &g);
    separator copy_spr(separator &source);
    void rclear_separator(separator &p);
    tree_property rfind_tree_properties(graph_type &g, int root);
    std::vector< std::vector <int> > rcreate_bfs_tree(graph_type &g, int root);
    int rfind_tallest_tree_root(graph_type &g);
    int find_farthest_node(graph_type &g, int v);
    std::vector<int> find_farthest_nodeall(graph_type &g, int v);
    solution rmake_solution(separator &spr);
    Graph_t make_boost_graph(graph_type &comp);
    graph_type make_graph_type(std::set <int> vertex_set);
};
