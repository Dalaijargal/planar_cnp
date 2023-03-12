#include "cndp_cycle_spr.h"
using namespace std;

using namespace std;
using namespace boost;

// struct for finding big face
struct face_visitor : public planar_face_traversal_visitor
{
	  void begin_face() {
	//	  std::cout << "New face: ";
		  tmp_face.resize(0);
		  counter = 0;
	  }
	  void end_face() {
	//	  std::cout<<"#: "<<counter;
		  if (counter > longest_face.size()){
			  longest_face = tmp_face;
		  }
	//	  std::cout << std::endl;
	  }

	  template <typename Vertex>
	  void next_vertex(Vertex v)
	  {
		  tmp_face.push_back(v);
	//	  std::cout << v << " ";
		  counter++;
	  }
	  std::vector<Vertex> get_lonsgest_face(){
		  return longest_face;
	  }
	  int counter;
	  std::vector<Vertex> longest_face, tmp_face;
};

// struct for finding tree edges
struct MyVisitor : boost::default_bfs_visitor {
    MyVisitor(EdgeSet& tree_edges) : tree_edges(tree_edges) {}

    void tree_edge(Edge e, const Graph_t& g) const {
        tree_edges.insert(e);
    }

    EdgeSet get_tree_edges(){
    	return tree_edges;
    }
  private:
    EdgeSet& tree_edges;
};

//constructor
cndp_cycle_spr::cndp_cycle_spr(int k1, int k2, bound_t bound, graph_type component, CNDP_Graph *graph_obj){
	this->total_k = k1;
	this->budget_k = k2;
	this->dim = graph_obj->dim;
	this->thecomponent = component;
	this->graph_object = graph_obj;
	this->bound = bound;
}

void cndp_cycle_spr::set_thecomponent(graph_type &c){
	this->thecomponent = c;
	//this->graph_object->graph_adj
}


vector<Vertex> cndp_cycle_spr::execute(){
		//vector <vector<Vertex> > cycle_spr_result;

		/*
		 * preparation part of fundamental cycle spr
		 */
		// convert graph_type graph to boost graph
	    Graph_t g1;
	    g1 = this->graph_object->make_boost_graph(this->thecomponent);

	    /*
	    std::vector<boost::graph_traits<Graph_t>::vertex_descriptor> vertices;
	    for (std::size_t i = 0; i < dim; ++i)
	      vertices.push_back(add_vertex(g1));
	    for (int i = 0; i < dim; i++){
	    	for (list<int>::iterator list_it = thecomponent.at(i).begin(); list_it != thecomponent.at(i).end(); list_it++)
	    		if (i < *list_it){
	    			boost::add_edge(vertices[i], vertices[*list_it], g1);
	    		}
	    }
	    */
//	    cout<<"num edges before tri: "<<num_edges(g1)<<endl;

	    // find biggest face
	    big_face_t big_face = find_biggest_face(g1);
	    Vertex root = big_face.root;
	    //root = height_minimize(g1, big_face);
	    //cout<<"Root in the primal: "<<root<<endl;

	    //cout<<"------big face--------"<<endl;
		//for (int i = 0; i < dim; i++){
		//	  if (big_face.table.at(i) == true){
		//		  cout<<i<<": "<<big_face.dist.at(i)<<endl;
		//	  }
		//} cout<<"\n";

	    // find tree-edges using BFS algorithm
	    std::set<Edge> tree_edges;
	    MyVisitor vis(tree_edges);
	    breadth_first_search(g1, root, boost::visitor(vis));
	    // get tree-edges and mark them in the graph
	    EdgeSet my_tree_edges = vis.get_tree_edges();
	    //cout<<"tree edges: "<<endl;
	    //for (EdgeSet::iterator it = my_tree_edges.begin(); it != my_tree_edges.end(); it ++){
	    //	cout<<*it<<endl;
	    //}

	    /*
	     * main part of computing fcs based on dual-tree and its properties
	     */

	    Graph_t g33;
	    g33 = triangulate(g1);
	    //cout<<"..."<<endl;
	    // redefine tree_edges in g33
	    set<Edge> my_tree_edges3;
	    my_tree_edges3 = redifine_tree_edges(my_tree_edges, g33);
	    //cout<<"...."<<endl;
	    // compute fcs
		vector <dual_tree_node> dual_tree;
	    dual_tree = compute_fcs(g33, my_tree_edges3, root);

	    //cout<<"....."<<endl;
		// select cycles that fit with my requirement;
	    vector <vector <Vertex> > selected_cycles = select_good_cycles(g33, dual_tree, this->bound, big_face, my_tree_edges3);

	    int len = INT_MAX;
	    int id = -1;
	    for (int i = 0; i < selected_cycles.size(); i ++){
	    	if (selected_cycles.at(i).size() < len){
	    		len = selected_cycles.at(i).size();
	    		id = i;
	    	}
	    }


	    /*
	     * examine the components divided by FC spr
	     */
	    /*
	    for (int i = 0; i < selected_cycles.size(); i++){
	    	vector<Vertex> thecycle = selected_cycles.at(i);
	    	vector<int> cycle_int(thecycle.size());
	    	for (int j = 0; j < thecycle.size(); j++){
	    		cycle_int.at(j) = thecycle.at(j);
	    		//cout<<cycle_int.at(j)<<" ";
	    	}
	    	separator myspr;
	    	myspr.vertex_list.push_back(cycle_int);
	    	solution sln = this->graph_object->rmake_solution(myspr);
	    	vector <int> index = this->graph_object->find_connected_components(sln.values);

	    	int num = *std::max_element(index.begin(),index.end());
	    	vector<int> counter(dim, 0);
	    	for (int t = 0; t < index.size(); t ++){
	    		//cout<<t<<": "<<index.at(t)<<endl;
	    		if (index.at(t) != -1)
	    			counter.at(index.at(t))++;
	    	}
	    	cout<<i<<"-th cycle's spr size: "<<cycle_int.size() <<"; # of components: "<<num+1<<":: "<<counter.at(0)<<"+"<<counter.at(1)<<"-"<<dim-counter.at(0)-counter.at(1)<<endl;
	    }
		*/
	    if (id != -1){
	    	//cout<<"id: "<<id<<"; # of cycles: "<<selected_cycles.size()<<endl;
	    	return selected_cycles.at(id);//cycle_spr_vec;
	    }
	    else{
	    	return vector <Vertex>();
	    	cout<<"Error in FCS"<<endl;
	    	exit(0);
	    }

} // end of execute

vector <vector <Vertex> > cndp_cycle_spr::select_good_cycles(Graph_t &g, vector <dual_tree_node> &mytree, bound_t bound, big_face_t bigFace, EdgeSet &tree_edges){
	/*
	 * preparation part: create tree-edge graph and tree-edge tree
	 */
	// create tree-edge boost graph
	Graph_t te_g;
	for (EdgeSet::iterator it = tree_edges.begin(); it != tree_edges.end(); it ++){
		Vertex trg = target(*it, g);
		Vertex src = source(*it, g);
		add_edge(src, trg, te_g);
	}

	// create a tree structure from tree_edges; each node store own parent in the primal
	// each node in te_tree stores only its parent and true if the node itself is on the big face, otherwise false
	int num = num_vertices(te_g);
	vector<pair<Vertex, bool> > te_tree;
	te_tree.resize(num);
    Vertex v_tmp;
    Vertex r = bigFace.root;
    vector<Vertex> rootLevel;  rootLevel.resize(0);
    vector<Vertex> tmpLevel;   tmpLevel.resize(0);
    vector<bool> visited(num, false);
    rootLevel.push_back(r);
    te_tree.at(r).first = r;
    te_tree.at(r).second = true;
    visited.at(r) = true;
    set<Vertex>::iterator set_it;
    typename boost::graph_traits<Graph_t>::out_edge_iterator oei, oei_end;
    // main loop for creating te_tree
    do{
    	for (vector<Vertex>::iterator it = rootLevel.begin(); it != rootLevel.end(); it ++){
//    		cout<<*it<<":-> "<<endl;
        	for(boost::tie(oei, oei_end) = boost::out_edges(*it, te_g); oei != oei_end; ++oei){
        		v_tmp = boost::target(*oei, te_g);
 //       		cout<<v_tmp<<endl;
        		if (visited.at(v_tmp) == false){
        			visited.at(v_tmp) = true;
 //       			cout<<v_tmp<<endl;
        			te_tree.at(v_tmp).first = *it;
        			te_tree.at(v_tmp).second = false;
        			if (bigFace.table.at(v_tmp)){
        				te_tree.at(v_tmp).second = true;
        			}
            		tmpLevel.push_back(v_tmp);
        		}
        	}
    	}
    	rootLevel.resize(0);
    	rootLevel = tmpLevel;
    	tmpLevel.resize(0);
    }while (rootLevel.size() != 0);
    //print out tree-edge tree data
    //for (int i = 0; i < te_tree.size(); i ++){
    //	cout<<"node id: "<<i<<" parent: "<<te_tree.at(i).first<<" on big face: "<<te_tree.at(i).second<<endl;
    //}
    //exit(0);



	/*
	 * main selection part that has three cases
	 */
	vector <vector <Vertex> > selected_cyc;
	//set<Vertex> tmp_cycle;
	set<Vertex> src_set;
	set<Vertex> trg_set;
	double up = 0.7;
	double lo = 0.3;
	int counter = 0;
	bool isSelected = false;
	//int diff = INT_MAX;
	//int balanced_cycle_id =-1;
	//cout<<"up: "<<dim*up<<"; low: "<<dim*lo<<endl;
	// Case 1: check if there is a good r2r cycle and restore it
	if (!isSelected){
		//cout<<"~~~~~~ offers for case 1: ~~~~~~~"<<endl;
		for (int i = 0 ; i < mytree.size(); i++){
			mytree.at(i).cycle_type = -1;
			int a = mytree.at(i).spr_info.out;
			int b = mytree.at(i).spr_info.in;
			int c = mytree.at(i).spr_info.spr;
			// basic criterion for case 1
			//up = 0.7;
			//lo = 0.3;
			if ((lo*dim < a) && (a < up*dim) && (lo*dim < b) && (b < up*dim)){
				// restore all cycles that meets basic requirement
				//Edge e1 = mytree.at(i).entrance_edge;
				//cout<<"entrance edge: "<<e1<<endl;
				set<Vertex> tmp_set;
				tmp_set = mytree.at(i).f_cycle;

				vector<Vertex> tmp_vec;
				int bf_counter = 0;
				for (set<Vertex>::iterator it = tmp_set.begin(); it != tmp_set.end(); it ++){
					//cout<<*it<<" ";
					if (bigFace.table.at(*it)){
						bf_counter ++;
					}
					tmp_vec.push_back(*it);
				}



				//cout<<endl;

				// check if the cycle is r2r.
				if (bf_counter < 2){
					//cout<<"bf_counter: "<<bf_counter<<" :=> ";
					mytree.at(i).cycle_type = 1;
					selected_cyc.push_back(tmp_vec);
					//for (set<Vertex>::iterator it = tmp_set.begin(); it != tmp_set.end(); it ++){
					//	cout<<*it<<" ";
					//}
					//cout<<endl;
					counter++;
				}
				//selected_cyc.push_back(tmp_set);
			}
		}
	}
	//cout<<"Case 1: selected cycle's number: "<<mytree.size()<<"-> "<<counter<<"\n\n\n";
	vector<int> vec;
	// Case 2: find beta-distanced cycles
	//if (!isSelected){
		//cout<<"\n\n~~~~~~ offers for case 2: ~~~~~~~"<<endl;
		for (int i = 0; i < mytree.size(); i++){
			int a = mytree.at(i).spr_info.out;
			int b = mytree.at(i).spr_info.in;
			int c = mytree.at(i).spr_info.spr;
			// basic criterion for case 1
			//up = 0.8;
			//lo = 0.2;
			// restore all cycles meeting the requirements
			if ((lo*dim < a) && (a < up*dim) && (lo*dim < b) && (b < up*dim) && (mytree.at(i).cycle_type = -1)){
				//cout<<"\n\nin: "<<mytree.at(i).spr_info.in<<"\n out: "<<mytree.at(i).spr_info.out<<"\n spr: "<<mytree.at(i).spr_info.spr<<endl;

				//cout<<"a b c: "<<a<<", "<<b<<", "<<c<<endl;
				set<Vertex> tmp_set = mytree.at(i).f_cycle;
				vector<Vertex> tmp_vec;
				// print out input cycle
				int bf_counter = 0;
				//cout<<"~~~~~~~~~ fcycle before cleaning: ";
				for (set<Vertex>::iterator it = tmp_set.begin(); it != tmp_set.end(); it ++){
					//cout<<*it<<" ";
					if (bigFace.table.at(*it)){
						bf_counter++;
					}
					tmp_vec.push_back(*it);
				}//cout<<endl;
				selected_cyc.push_back(tmp_vec);
				//cout<<"big face size: "<<bf_counter<<endl;
/*
				// check if the cycle is not r2r one?
				if (bf_counter >= 2){
					//cout<<"# of big face nodes: "<<bf_counter<<": "<<"path spr"<<endl;
					//for (set<Vertex>::iterator it = tmp_set.begin(); it != tmp_set.end(); it ++){
					//	cout<<*it<<" ";
					//} cout<<endl;

					vector<Vertex> path1;
					vector<Vertex> path2;

					// cleaning this spr using tree-edge graph and big face table and dist
					Edge e1 = mytree.at(i).entrance_edge;
					//cout<<"entrance edge: "<<e1<<endl;
					Vertex src = source(e1, g);
					Vertex tmp_src = src;
					Vertex trg = target(e1, g);

					// examine first path starting from src of enterance edge
					bool onPath = false;
					do{
						Vertex parent = te_tree.at(tmp_src).first;
						if ((!bigFace.table.at(parent)) && (bigFace.table.at(tmp_src))){
							onPath = true;
						}
						if ((bigFace.table.at(parent)) && (!bigFace.table.at(tmp_src)) && (onPath)){
							onPath = false;
							//cout<<tmp_src<<" "<<parent;
							path1.push_back(tmp_src);
							path1.push_back(parent);
							break;
						}
						if (onPath){
						//	cout<<tmp_src<<" ";
							path1.push_back(tmp_src);
						}
						set<Vertex>::iterator set_it = tmp_set.find(parent);
						if ((set_it == tmp_set.end())||(tmp_src == parent))
							break;
						tmp_src = parent;
					}while(1);

					// examine second path starting from target of entrance edge
					tmp_src = trg;
					onPath = false;
					do{
						Vertex parent = te_tree.at(tmp_src).first;
						if ((!bigFace.table.at(parent)) && (bigFace.table.at(tmp_src))){
							onPath = true;
						}
						// break for one path spr
						if ((bigFace.table.at(parent)) && (!bigFace.table.at(tmp_src)) && (onPath)){
							onPath = false;
							//cout<<tmp_src<<" "<<parent;
							path2.push_back(tmp_src);
							path2.push_back(parent);
							break;
						}
						if (onPath){
							//cout<<tmp_src<<" ";
							path2.push_back(tmp_src);
						}
						// break for two path spr
						set<Vertex>::iterator set_it = tmp_set.find(parent);
						if ((set_it == tmp_set.end())||(tmp_src == parent)){
							break;
						}
						//if (tmp_src == parent)
						//	break;

						tmp_src = parent;
					}while(1);
			//exit(0);
					// chose one from two path
					if ((path1.size() != 0) && (path2.size() != 0)){
						path2.erase(--path2.end());//path2.end();

						//cout<<"both path has discovered spr"<<endl;
						//cout<<"spr discovered from path1: ";
						//for (int j = 0; j < path1.size(); j++)
						//	cout<<path1.at(j)+1<<" ";

						//cout<<endl<<"spr discovered from path2: ";
						//for (int j = 0; j < path2.size(); j++)
						//	cout<<path2.at(j)+1<<" ";
						//cout<<endl;

						for (int j = path2.size() - 1; j >= 0; j --){
							path1.push_back(path2.at(j));
						}

						//cout<<endl<<"merged path: ";
						//for (int j = 0; j < path1.size(); j++)
						//	cout<<path1.at(j)+1<<" ";
						//cout<<endl;
						selected_cyc.push_back(path1);
					}
					else{
						if (path1.size() != 0){
							// examine distance from root, choose it if it meets the requirement

							selected_cyc.push_back(path1);
							//cout<<"spr discovered from path1: ";
							//for (int j = 0; j < path1.size(); j++)
							//	cout<<path1.at(j)<<" ";
							//cout<<"; end points: "<<*path1.begin()<<" <-> "<<*--path1.end()<<endl;
							//for ()
						}
						else if (path2.size() != 0){
							// examine distance from root, choose it if it meets the requirement

							selected_cyc.push_back(path2);
							//cout<<"spr discovered from path2: ";
							//for (int j = 0; j < path2.size(); j++)
							//	cout<<path2.at(j)<<" ";
							//cout<<"; end points: "<<*path2.begin()<<" <-> "<<*--path2.end()<<endl;
						}
					}

					//cout<<endl<<"---------"<<endl;
				}// end of path spr cleaning and restoring
*/

			}// end of basic criteria
		}// end of for

		if (selected_cyc.size() > 0)
			isSelected = true;
	//}// end of case 2

	// Case 3: return normal spr
	//if (!isSelected){
	//	cout<<"\n\n~~~~~~ offers for case 3: ~~~~~~~"<<endl;
	//}


	return selected_cyc;
}

int cndp_cycle_spr::find_height(Graph_t &g, Vertex v){

	//cout<<"hello from find height"<<endl;
    vector<Vertex> rootLevel;  rootLevel.resize(0);
    vector<Vertex> tmpLevel;   tmpLevel.resize(0);
    rootLevel.push_back(v);
    vector<bool> visited(dim, false);
    Vertex v_tmp;
    rootLevel.push_back(v);
    visited.at(v) = true;
    int height = 0;

    //set<Vertex>::iterator set_it;
    typename boost::graph_traits<Graph_t>::out_edge_iterator oei, oei_end;
    do{
    	for (vector<Vertex>::iterator it = rootLevel.begin(); it != rootLevel.end(); it ++){
    		//cout<<*it<<":-> "<<endl;
    		//visited.at(*it) = true;
        	for(boost::tie(oei, oei_end) = boost::out_edges(*it, g); oei != oei_end; ++oei)
        	{
        		v_tmp = boost::target(*oei, g);
 //       		cout<<v_tmp<<endl;
        		if (visited.at(v_tmp) == false){
        			//cout<<v_tmp<<endl;
            		tmpLevel.push_back(v_tmp);
            		visited.at(v_tmp) = true;
        		}
        	}
        	//cout<<endl;
    	}
    	rootLevel.resize(0);
    	rootLevel = tmpLevel;
    	tmpLevel.resize(0);
    	height ++;
    }while (rootLevel.size() != 0);


	return height;
} // end of find_height

// find a root on big face that cause a short bfs tree
Vertex cndp_cycle_spr::height_minimize(Graph_t &g, set<Vertex> &big_face){
	Vertex root;
	int r = rand() % big_face.size();
	int min_height = INT_MAX;
	Vertex short_root;
	int counter = 0;
	for (set<Vertex>::iterator it = big_face.begin(); it != big_face.end(); it ++){
		root = *it;
		int h = find_height(g, root);
		if (min_height > h){
			min_height = h;
			short_root = root;
		}
		counter++;
//		if (counter > 2)
//			break;
//		cout<<"root: "<<root<<" "<<h<<endl;
	}
//	cout<<"selected root: "<<short_root<<" "<<min_height<<endl;
	return short_root;
}

// find big face
// it is slow component
big_face_t cndp_cycle_spr::find_biggest_face(Graph_t &g_orig){
	  // Initialize the interior edge index
	  Graph_t g;
	  g = g_orig;
	  /*
	   * find vertices on big face
	   */
	  property_map<Graph_t, edge_index_t>::type e_index = get(edge_index, g);
	  graph_traits<Graph_t>::edges_size_type edge_count = 0;
	  graph_traits<Graph_t>::edge_iterator ei, ei_end;
	  for(tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
	    put(e_index, *ei, edge_count++);
	  // Test for planarity - we know it is planar, we just want to
	  // compute the planar embedding as a side-effect
	  typedef std::vector< graph_traits<Graph_t>::edge_descriptor > vec_t;
	  std::vector<vec_t> embedding(num_vertices(g));
	  if (boyer_myrvold_planarity_test(boyer_myrvold_params::graph = g,
	                                   boyer_myrvold_params::embedding =
	                                       &embedding[0]
	                                   )
	      )
	    std::cout << "";// "Input graph is planar" << std::endl;
	  else{
	    std::cout << "Input graph is not planar" << std::endl;
	    exit(0);
	  }
	  face_visitor v_vis;
	  planar_face_traversal(g, &embedding[0], v_vis);
	  std::set<Vertex> big_face_vertex;
	  std::vector<Vertex> bigt_face_v = v_vis.get_lonsgest_face();
	  //cout<<"longest face's length: "<<bigt_face.size()<<endl;
	  for (int i = 0; i < bigt_face_v.size(); i++){
		  big_face_vertex.insert(bigt_face_v.at(i));
	  }


	  // find root from big face
	  Vertex root = height_minimize(g_orig, big_face_vertex);


	  //?????
	  //root = vertex(93, g);
	  //cout<<"root: "<<root<<endl;

	  /*
	   * create big face table and distance from root
	   */
	  //pair<vector <bool>, vector <int> > bf_data;
	  big_face_t bf_dat;
	  bf_dat.length = big_face_vertex.size();
	  // find big face table
	  //bf_data.first.resize(dim, false);
	  bf_dat.table.resize(dim, false);
	  for (set<Vertex>::iterator it = big_face_vertex.begin(); it != big_face_vertex.end(); it ++){
		  //bf_data.first.at(*it) = true;
		  bf_dat.table.at(*it) = true;
	  }

	  // create big face tree

	  Vertex v_tmp;
	  vector<Vertex> rootLevel;  rootLevel.resize(0);
	  vector<Vertex> tmpLevel;   tmpLevel.resize(0);
	  vector<bool> visited(dim, false);
	  rootLevel.push_back(root);
	  visited.at(root) = true;
	  bf_dat.root = root;
	    //int level = 1;
	  set<Vertex>::iterator set_it;
	  typename boost::graph_traits<Graph_t>::out_edge_iterator oei, oei_end;

	  // temporary remove one of child of root
	  //cout<<"nodes on big face before removing: "<<endl;
	  //for (set<Vertex>::iterator it = big_face.begin(); it != big_face.end(); it ++){
	//	  cout<<*it<<" "<<bf_dat.table.at(*it)<<endl;
	  //}	  cout<<endl;
	  Vertex v_removed;
	  for(boost::tie(oei, oei_end) = boost::out_edges(root, g_orig); oei != oei_end; ++oei)  {
	  		v_removed = boost::target(*oei, g_orig);
	  		if (bf_dat.table.at(v_removed)) {
	  			bf_dat.table.at(v_removed) = false;
	  			//big_face.erase(v_removed);
	  			//cout<<"removed v: "<<v_removed<<endl;
	  			break;
	  		}
	  }

//	  cout<<"nodes on big face after removing: "<<endl;
//	  for (set<Vertex>::iterator it = big_face.begin(); it != big_face.end(); it ++){
//		  cout<<*it<<" "<<bf_dat.table.at(*it)<<endl;
//	  }	  cout<<endl;

	  //cout<<"big face tree series: "<<endl;
	  int dist = 0;
	  //bf_data.second.resize(dim, -1);
	  //bf_data.second.at(root) = dist;
	  bf_dat.dist.resize(dim, -1);
	  bf_dat.dist.at(root) = dist;
	  do{
		  for (vector<Vertex>::iterator it = rootLevel.begin(); it != rootLevel.end(); it ++){
	    		//cout<<*it<<":-> "<<endl;
	        	for(boost::tie(oei, oei_end) = boost::out_edges(*it, g_orig); oei != oei_end; ++oei)
	        	{
	        		//cout<<*oei<<" ";
	        		v_tmp = boost::target(*oei, g_orig);
	        		//cout<<v_tmp<<endl;
	        		//set<Vertex>::iterator it_set;
	        		//it_set = big_face.find(v_tmp);
	        		// examine if it is on big face
	        		//if (it_set != big_face.end()){

	        		if (bf_dat.table.at(v_tmp)){
		        		if (visited.at(v_tmp) == false){
		        			//exit(0);
		        			dist++;
		        			visited.at(v_tmp) = true;
		        			//cout<<dist<<": "<<v_tmp<<endl;
		        			bf_dat.dist.at(v_tmp) = dist;
		        			//bf_data.second.at(v_tmp) = dist;

		            		tmpLevel.push_back(v_tmp);
		        		}
	        		}

	        	}
	    	}
	    	rootLevel.resize(0);
	    	rootLevel = tmpLevel;
	    	tmpLevel.resize(0);
	  }while (rootLevel.size() != 0);
	 //exit(0);
	  //bf_data.second.at(v_removed) = ++dist;
	  bf_dat.table.at(v_removed) = true;
	  bf_dat.dist.at(v_removed) = ++dist;
	  //for (int i = 0; i < dim; i++){
		//  if (bf_dat.table.at(i) == true){
		//	  cout<<i<<": "<<bf_dat.dist.at(i)<<endl;
		 // }
	  //}

return bf_dat;

}

// find common vertex when 2 cycle merged
Vertex cndp_cycle_spr::find_o(vector<pair<Vertex, bool> > &te_tree, set<Vertex> &fc_intersect, Vertex p2, Vertex &root_g){
	//cout<<"find o: "<<p2<<"-> ";
	for (int i = 0; i < te_tree.size(); i++){
		te_tree.at(i).second = false;
	}

	set<Vertex>::iterator set_it;
	for (set_it = fc_intersect.begin(); set_it != fc_intersect.end(); set_it++){
		te_tree.at(*set_it).second = true;
	}

	Vertex p22 = p2;
	Vertex parent = p2;
	bool status = false;
	int counter = 0;
	do{
//		cout<<p2<<"-> "<<parent<<"= parent status: "<<status<<endl;
		p2 = parent;
		parent = te_tree.at(p2).first;
		status = te_tree.at(parent).second;
		counter++;
	}while ( status && (parent != root_g));

//	cout<<"~~~distance from "<<p22<<" to "<<parent<<": "<<counter<<endl;
//	exit(0);
	if (!status){
		//cout<<p2<<endl;
		return p2;
	}
	else{
		//cout<<parent<<endl;
		return parent;
	}

}

vector<dual_tree_node> cndp_cycle_spr::compute_fcs(Graph_t &g, EdgeSet &tree_edges, Vertex root_g){

	/*
	 * find dual graph tree and its properties
	 */
    // mark tree and non-tree edges on g
    gt_edgeIterator ei, ei_end;
    pm_edgePropertyType pm_edgeWeight = get(edge_weight, g);
    for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei){
    	put(pm_edgeWeight, *ei, 0);
    }
    std::set<Edge>::iterator edgeIterator;
    for (edgeIterator = tree_edges.begin(); edgeIterator != tree_edges.end(); edgeIterator++){
    	put(pm_edgeWeight, *edgeIterator, 1);
    }

    /*
     *  create dual graph and rank non-tree edges
     */
    vector <dual_tree_node> tree;
    tree = rank_ntree_edges(root_g, g, tree_edges);

	// compute vertex id map
	vector <int> id_map;
	id_map.resize(tree.size());
	for (int i = 0; i < tree.size(); i++){
		id_map[tree.at(i).vertex_id] = i;
	}

	/*
	 * create a graph with tree_edges in the primal to calculate op2
	 */

	Graph_t te_g;
//	int num = num_vertices(te_g);
	for (EdgeSet::iterator it = tree_edges.begin(); it != tree_edges.end(); it ++){
		Vertex trg = target(*it, g);
		Vertex src = source(*it, g);
		add_edge(src, trg, te_g);
	}

	// create a tree structure from tree_edges; each node store own parent in the primal
	// each node in te_tree stores only its parent and true if the parent is in the cycle
	int num = num_vertices(te_g);
	vector<pair<Vertex, bool> > te_tree;
	te_tree.resize(num);
    Vertex v_tmp;
    vector<Vertex> rootLevel;  rootLevel.resize(0);
    vector<Vertex> tmpLevel;   tmpLevel.resize(0);
    vector<bool> visited(num, false);
    rootLevel.push_back(root_g);
    visited.at(root_g) = true;
    //int level = 1;
    set<Vertex>::iterator set_it;
    typename boost::graph_traits<Graph_t>::out_edge_iterator oei, oei_end;
    do{
    	for (vector<Vertex>::iterator it = rootLevel.begin(); it != rootLevel.end(); it ++){
//    		cout<<*it<<":-> "<<endl;
        	for(boost::tie(oei, oei_end) = boost::out_edges(*it, te_g); oei != oei_end; ++oei)
        	{
        		v_tmp = boost::target(*oei, te_g);
 //       		cout<<v_tmp<<endl;
        		if (visited.at(v_tmp) == false){
        			visited.at(v_tmp) = true;
 //       			cout<<v_tmp<<endl;
        			te_tree.at(v_tmp).first = *it;
            		tmpLevel.push_back(v_tmp);
        		}
        	}
    	}
    	rootLevel.resize(0);
    	rootLevel = tmpLevel;
    	tmpLevel.resize(0);
    }while (rootLevel.size() != 0);

    /*
     * Main part of computing fcs
     */

	int N = num_vertices(g);
	for (int i = 0; i < tree.size() - 1; i++){
		//cout<<"-----------------------"<<endl;
		//cout<<"entrance edge: "<<tree.at(i).entrance_edge<<endl;
		// if the node is leaf, ...
		if (tree.at(i).child_id.size() == 0){
			tree.at(i).spr_info.in = 0;
			tree.at(i).spr_info.spr = 3;
			tree.at(i).spr_info.out = N - 3 - 0;
			for (int j = 0; j < 3; j++){
				tree.at(i).f_cycle.insert( tree.at(i).face_vertex_list.at(j) );
//				cout<<tree.at(i).face_vertex_list.at(j)<<endl;
			}
			//cout<<"hello from case leaf: "<< tree.at(i).entrance_edge <<endl;
		}
		// case A
		else if (tree.at(i).child_id.size() == 1){
			int a, b, c;
			int child_id = *tree.at(i).child_id.begin();
			a = tree.at(id_map[child_id]).spr_info.out;
			b = tree.at(id_map[child_id]).spr_info.in;
			c = tree.at(id_map[child_id]).spr_info.spr;

			Edge e = tree.at(i).entrance_edge;
			Vertex p2, p3;
			p2 = source(e, g);
			p3 = target(e, g);

			set<Vertex>::iterator vi1, vi2;
			set<Vertex> childs_cycle;

			childs_cycle = tree.at(id_map[child_id]).f_cycle;
			vi1 = childs_cycle.find(p2);
			vi2 = childs_cycle.find(p3);

			// case A1
			if ( ( (vi1 != childs_cycle.end()) && (vi2 == childs_cycle.end()) )
			|| ( (vi1 == childs_cycle.end()) && (vi2 != childs_cycle.end()) ) ){
				tree.at(i).spr_info.out = a - 1;
				tree.at(i).spr_info.in  = b;
				tree.at(i).spr_info.spr = c + 1;
				// assign vertices in fundamental cycle
				for (set<Vertex>::iterator it = tree.at(id_map[child_id]).f_cycle.begin();
						it != tree.at(id_map[child_id]).f_cycle.end(); it ++){
					tree.at(i).f_cycle.insert(*it);
				}
				tree.at(i).f_cycle.insert(p2);
				tree.at(i).f_cycle.insert(p3);
				//cout<<"hello from case A1: "<<e<<endl;
			}

			// case A2
			if ((vi1 != childs_cycle.end()) && (vi2 != childs_cycle.end())) {
				tree.at(i).spr_info.out = a;
				tree.at(i).spr_info.in = b + 1;
				tree.at(i).spr_info.spr = c - 1;
				// assign vertices in fundamental cycle
				for (set<Vertex>::iterator it = tree.at(id_map[child_id]).f_cycle.begin();
						it != tree.at(id_map[child_id]).f_cycle.end(); it ++){
					tree.at(i).f_cycle.insert(*it);
				}
				for (int j = 0; j < 3; j++){
					Vertex v = tree.at(i).face_vertex_list.at(j);
					if ((p2 != v)&&(p3 != v)){
						tree.at(i).f_cycle.erase(v);
						break;
					}
				}
				//cout<<"hello from case A2: "<<e<<endl;
			}
		}
		// case B
		else if (tree.at(i).child_id.size() == 2){
			int a1, a2, b1, b2, c1, c2;
			int child_id1 = *tree.at(i).child_id.begin();
			int child_id2 = *++tree.at(i).child_id.begin();
//			cout<< id_map[child_id1] << " " << id_map[child_id2]<<endl;

			a1 = tree.at(id_map[child_id1]).spr_info.out;
			b1 = tree.at(id_map[child_id1]).spr_info.in;
			c1 = tree.at(id_map[child_id1]).spr_info.spr;

			a2 = tree.at(id_map[child_id2]).spr_info.out;
			b2 = tree.at(id_map[child_id2]).spr_info.in;
			c2 = tree.at(id_map[child_id2]).spr_info.spr;

			Vertex p2;// = vertex(4, g);
			Edge e1 = *tree.at(i).outgoing_edges.begin();
			Edge e2 = *++tree.at(i).outgoing_edges.begin();
			//cout<<"outgoing edges: "<<e1<<" "<<e2<<endl;
			Vertex src1, src2, trg1, trg2;
			src1 = source(e1, g);
			trg1 = target(e1, g);
			src2 = source(e2, g);
			trg2 = target(e2, g);
			if (src1 == src2)
				p2 = src1;
			else if (trg1 == trg2)
				p2 = trg1;
			else if (src1 == trg2)
				p2 = trg2;
			else if (src2 == trg1)
				p2 = trg1;

			// assign vertices in fundamental cycle
			set<Vertex> fc1 = tree.at(id_map[child_id1]).f_cycle;
			set<Vertex> fc2 = tree.at(id_map[child_id2]).f_cycle;
			set<Vertex> fc_intersection;
			std::set_intersection(fc1.begin(), fc1.end(),
								fc2.begin(), fc2.end(),
								std::inserter(fc_intersection, fc_intersection.begin()));
			set<Vertex> fc_union;

			for (set<Vertex>::iterator it = fc1.begin(); it != fc1.end(); it ++){
				fc_union.insert(*it);
			}
			for (set<Vertex>::iterator it = fc2.begin(); it != fc2.end(); it ++){
				fc_union.insert(*it);
			}
			Vertex o = find_o(te_tree, fc_intersection, p2, root_g);
		    std::set<Vertex> fc_diff;

		    std::set_difference(fc_union.begin(), fc_union.end(), fc_intersection.begin(), fc_intersection.end(),
		                        std::inserter(fc_diff, fc_diff.begin()));

		    for (set<Vertex>::iterator it = fc_diff.begin(); it != fc_diff.end(); it ++)
		    	tree.at(i).f_cycle.insert(*it);
		    tree.at(i).f_cycle.insert(o);

			int op2 = fc_intersection.size(); //dist_te[p2];
			tree.at(i).spr_info.in = b1 + b2 + op2 -1;
			tree.at(i).spr_info.spr = (c1 - op2) + (c2 - op2) + 1;
			tree.at(i).spr_info.out = N - (tree.at(i).spr_info.in+tree.at(i).spr_info.spr);

			//cout<<"hello from case B: "<<tree.at(i).entrance_edge<<endl;
			//cout<<"child cycles: ";
			//for (set<Vertex>::iterator it = fc1.begin(); it != fc1.end(); it ++){
			//	cout<<*it<<" ";
			//} cout<<"| ";
			//for (set<Vertex>::iterator it = fc2.begin(); it != fc2.end(); it ++){
			//	cout<<*it<<" ";
			//} cout<<endl;

		}// end case B
		//cout<<"cycle: ";
		//for (set<Vertex>::iterator it = tree.at(i).f_cycle.begin(); it != tree.at(i).f_cycle.end(); it ++){
		//	cout<<*it<<" ";
		//} cout<<endl;

	}

return tree;
} // end compute fcs

bool wayToSort(dual_tree_node a, dual_tree_node b) {
	return a.bfs_level > b.bfs_level;
}

vector <dual_tree_node> cndp_cycle_spr::rank_ntree_edges(Vertex root, Graph_t &g, EdgeSet &my_tree){
    /*
     * find edges of dual graph on non-tree-edges: it will be the dual tree
     */

	// Create an empty graph to hold the dual
    Graph_t dual_g;
    //Initialize the interior edge index
    pm_edgeIndexType e_index = get(edge_index, g);
    graph_traits<Graph_t>::edges_size_type edge_count = 0;
    graph_traits<Graph_t>::edge_iterator ei, ei_end;
    for(tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
    {
      put(e_index, *ei, edge_count++);
    }

    // Compute the planar embedding - we know the input graph is planar,
    // so we're ignoring the return value of the test
    typedef std::vector< graph_traits<Graph_t>::edge_descriptor > vec_t;
    std::vector<vec_t> embedding(num_vertices(g));
    boyer_myrvold_planarity_test(boyer_myrvold_params::graph = g,
                                 boyer_myrvold_params::embedding = &embedding[0]
                                 );

    typedef typename std::vector <std::vector <Vertex> > face_vector_t;
    face_vector_t::iterator fvIter;
    std::vector <Vertex>::iterator vertexIter;

    // create dual graph from g
    face_vector_t faceVector = create_dual_graph(g, dual_g, &embedding[0], get(edge_index,g), get(edge_weight,g));

    //cout<<"edges in dual graph"<<endl;
    //for (tie(ei, ei_end) = edges(dual_g); ei != ei_end; ei++){
    //	cout<<*ei<<endl;
    //}
    //cout<<"faces and id in dual g"<<endl;
    //for (int i = 0; i < faceVector.size(); i ++){
    //	cout<<i<<": ";
    //	for (int j = 0; j < faceVector.at(i).size(); j++){
    //		cout<<faceVector.at(i).at(j)<<" ";
    //	}
    //	cout<<endl;
    //}

    /*
     * make a face including root the root face in dual graph
     */
	//define root of dual graph tree
	    int root_face_id = 0;
	    for (fvIter = faceVector.begin(); fvIter != faceVector.end(); fvIter++){
	  	  if(std::find(fvIter->begin(), fvIter->end(), root) != fvIter->end()) {
	  		  //std::cout<<"Root face: "<<fvIter->at(0)<<" "<<fvIter->at(1)<<" "<<fvIter->at(2)<<" "<<std::endl;
	  		  break;
	  	  }
	  	  root_face_id++;
	    }
	    //cout<<"root face id: "<<root_face_id<<endl;

	/*
	 * construct dual tree; vertex id, face_vertex_list, face_edge_list,
	 */

	    int num = num_vertices(dual_g);
	    //cout<<"# of vertices in dual: "<<num<<endl;
	    vector <dual_tree_node> dual_tree(num);
	    for (int i = 0; i < num; i++){
	    	//	  	  cout<<i<<endl;
	    		  	  // set dual_tree.vertex_id
	    		  	  dual_tree.at(i).vertex_id = vertex(i, dual_g);
	    	//	  	  cout<<faceVector.at(i).at(0)<<" "<<faceVector.at(i).at(1)<<" "<<faceVector.at(i).at(2)<<" "<<endl;

	    		  	  // set dual_tree.face_edge_list and its status; and face_vertex_list
	    		  	  std::set<Edge>::iterator eIt;
	    		  	  Vertex v0 = faceVector.at(i).at(0);
	    		  	  Vertex v1 = faceVector.at(i).at(1);
	    		  	  Vertex v2 = faceVector.at(i).at(2);
	    		  	  dual_tree.at(i).face_vertex_list.push_back(v0);
	    		  	  dual_tree.at(i).face_vertex_list.push_back(v1);
	    		  	  dual_tree.at(i).face_vertex_list.push_back(v2);
	    		  	  bool f1 = false, f2 = false, f3 = false;


	    		  	  // define face edges with face_edge_list
	    		  	  Edge e1, e2;
	    		  	  dual_tree.at(i).face_edge_list.resize(3);
	    	//	  	  set<Edge>::iterator eIt2;
	    		  	  EdgeSet::iterator eIt2;
	    		  	  EdgeSet::iterator eIt1;
	    		  	  boost::tie(e1, f1) = edge(v0, v1, g);
	    		  	  boost::tie(e2, f1) = edge(v1, v0, g);
	    		  	  if (f1){
	    		//  		cout<<"-"<<e1<<endl;
	    		  		  dual_tree.at(i).face_edge_list.at(0).edge.first =e1;
	    		  		  eIt1 = my_tree.find(e1);
	    		  		  eIt2 = my_tree.find(e2);

	    		  		  if ((eIt1 != my_tree.end())||(eIt2 != my_tree.end()))
	    		  			  dual_tree.at(i).face_edge_list.at(0).edge.second = -1;
	    		  		  else
	    		  			  dual_tree.at(i).face_edge_list.at(0).edge.second = 1;
	    		  	  }
	    		  	  else{
	    		  		  cout<<"There is not that edge"<<endl;
	    		  		  exit(0);
	    		  	  }

	    		  	  boost::tie(e1, f2) = edge(v1, v2, g);
	    		  	  boost::tie(e2, f2) = edge(v2, v1, g);
	    		  	  if (f2){
	    		  		  dual_tree.at(i).face_edge_list.at(1).edge.first = e1;
	    		  		  eIt1 = my_tree.find(e1);
		    		  	  eIt2 = my_tree.find(e2);
		    		  	  if ( (eIt1 != my_tree.end()) || (eIt2 != my_tree.end()) )
		    				dual_tree.at(i).face_edge_list.at(1).edge.second = -1;
		    	 		  else
		    	  			dual_tree.at(i).face_edge_list.at(1).edge.second = 1;
	    		  	  }
	    		  	  else{
	    		  		  cout<<"There is not that edge"<<endl;
	    		  		  exit(0);
	    		  	  }

	    		  	  boost::tie(e1, f3) = edge(v2, v0, g);
	    		  	  boost::tie(e2, f3) = edge(v0, v2, g);

	    		  	  if (f3){
	    			  		dual_tree.at(i).face_edge_list.at(2).edge.first = e1;
	    			  		eIt1 = my_tree.find(e1);
	    			  		eIt2 = my_tree.find(e2);
	    			  		if ( (eIt1 != my_tree.end()) || (eIt2 != my_tree.end()) )
	    			  			dual_tree.at(i).face_edge_list.at(2).edge.second = -1;
	    			  		else
	    			  			dual_tree.at(i).face_edge_list.at(2).edge.second = 1;

	    		  	  }
	    		  	  else{
	    		  		  cout<<"There is not that edge"<<endl;
	    		  		  exit(0);
	    		  	  }

	    		  	  // assign bfs level
	    		  	  //set entrance-edge

	    }// end construction dual-tree


		/*
		 * construct dual tree, continue; bfs_level, parent_id, child_id, entrance_edge, outgoing edge
		 * ***/

//	      cout<<"------ bfs level ------"<<endl;
	      typename boost::graph_traits<Graph_t>::out_edge_iterator oei, oei_end;
	      typename boost::graph_traits<Graph_t>::vertex_iterator vi, vi_end;
	      typename boost::graph_traits<Graph_t>::vertex_descriptor v1, v2;
	      v1 = vertex(0, dual_g);

	      Vertex R = vertex(root_face_id, dual_g);
	      Vertex v_tmp;
	      vector <int> bfs_level;
	      int level = 0;
	      vector<Vertex> rootLevel;
	      rootLevel.resize(0);
	      vector<Vertex> tmpLevel;
	      tmpLevel.resize(0);
	      vector<Vertex>::iterator it;
	      vector<bool> visited(num, false);
	      rootLevel.push_back(R);
	      visited.at(R) = true;
//	      cout<<"Root: "<<level<<": "<<R<<endl;
	      // traverse bfs tree
	      do{
		      	level++;
		      	for (it = rootLevel.begin(); it != rootLevel.end(); it ++){
		  //    		cout<<*it<<": ";


		          	for(boost::tie(oei, oei_end) = boost::out_edges(*it, dual_g); oei != oei_end; ++oei)
		          	{
		          		v_tmp = boost::target(*oei, dual_g);
		          		if (visited.at(v_tmp)==false){
		          			visited.at(v_tmp) = true;
		              		dual_tree.at(v_tmp).parent_id = *it;
		              		dual_tree.at(*it).child_id.push_back(v_tmp);
//		              		cout<<level<<": "<<*it<<"-> "<<v_tmp<<endl;
		              		dual_tree.at(v_tmp).bfs_level = level;
		              		tmpLevel.push_back(v_tmp);
		              		// compute entrance_edge; and outgoing_edge
		              		// every connected 2 face has only one non-tree-edge
		              		edge_t e1, e2;
		              		bool found = false;
		              		// iterate all edges in outgoing face
		              		for (int i = 0; i < 3; i++){
		              			e1 = dual_tree.at(*it).face_edge_list.at(i);
		              			if (e1.edge.second == 1){
//		              				cout<<e1.edge.first<<"-> "<<endl;
		   				    		Vertex src1, trg1, src2, trg2;
		   				    		src1 = source(e1.edge.first, g);
		   				    		trg1 = target(e1.edge.first, g);

		   				    		// iterate all edges in the entering face
		              				for (int j = 0; j < 3; j++){
		              					e2 = dual_tree.at(v_tmp).face_edge_list.at(j);
		              				     if (e2.edge.second == 1){
		              				    	//cout<<e2.edge.first<<" "<<endl;
		           				    		src2 = source(e2.edge.first, g);
		           				    		trg2 = target(e2.edge.first, g);
		           				    		if (((src1 == src2)&&(trg1 == trg2))||((src1 == trg2)&&(trg1 == src2))){
		           				    			found = true;
//		           				    			cout<<"entrance edge: "<<dual_tree.at(v_tmp).face_edge_list.at(j).edge.first<<endl;
		           				    			dual_tree.at(v_tmp).entrance_edge = e2.edge.first;
		           				    			dual_tree.at(*it).face_edge_list.at(i).edge.second = 2;
		           				    			dual_tree.at(v_tmp).face_edge_list.at(j).edge.second = 2;
		           				    		}
		           				    		//else{
	//	           				    		//	cout<<"outgoing edge: "<<e2.edge.first<<" ";
		           				    		//	dual_tree.at(v_tmp).outgoing_edges.push_back(e2.edge.first);
		           				    		//}
	//	           				    		cout<<endl;
		              				     }
				              		}
		              			}
		              			if (found){
		   				    		// iterate all edges in the entering face to mark its outgoing edges
		              				// each face has only one entrance edge and so other non-tree edges are outgoing
		              				for (int j = 0; j < 3; j++){
		              					e2 = dual_tree.at(v_tmp).face_edge_list.at(j);
		              				     if (e2.edge.second == 1){
		           				    		dual_tree.at(v_tmp).outgoing_edges.push_back(e2.edge.first);
		              				     }
				              		}

		              				break;
		              			}
		              		}

		          		}
		          	}

		      	}
		      	//cout<<endl;
		      	rootLevel.resize(0);
		      	rootLevel = tmpLevel;
		      	tmpLevel.resize(0);

	      }while (rootLevel.size() != 0); // end of computing bfs level

	      // sort dual_tree in terms to bfs_level
	      sort(dual_tree.begin(), dual_tree.end(), wayToSort);

//	      cout<<"after sorting in terms to bfs_level"<<endl;
	      //print_dual_tree(dual_tree);
	      //exit(0);

return dual_tree;
}

// re-define tree edges after triangulation because the graph has changed
set<Edge> cndp_cycle_spr::redifine_tree_edges(set<Edge> &tree_edges, Graph_t &g3){
    Vertex v1, v2;
    set<Edge> my_tree_edges3;
    Edge e;
    bool b;
    for (std::set<Edge>::iterator edgeIterator = tree_edges.begin(); edgeIterator != tree_edges.end(); edgeIterator++){
    	v1 = source(*edgeIterator, g3);
    	v2 = target(*edgeIterator, g3);
    	boost::tie(e, b) = edge(v1, v2, g3);
    	if (b){
    		my_tree_edges3.insert(e);
    	}
    	else {
    		cout<<"there is not that edge in g3"<<endl;
    		exit(0);
    	}
    }

    // mark tree and non-tree edges on g3, not on g
    gt_edgeIterator ei, ei_end;
    pm_edgePropertyType pm_edgeWeight = get(edge_weight, g3);
    for (tie(ei, ei_end) = edges(g3); ei != ei_end; ++ei){
    	put(pm_edgeWeight, *ei, 0);
    }
    for (std::set<Edge>::iterator edgeIterator = my_tree_edges3.begin(); edgeIterator != my_tree_edges3.end(); edgeIterator++){
    	put(pm_edgeWeight, *edgeIterator, 1);
    }
    return my_tree_edges3;
}


template <typename Graph_t>
void cndp_cycle_spr::print_graph(Graph_t& g)
{
  std::cout << "vertices: " << std::endl;
  typename graph_traits<Graph_t>::vertex_iterator vi, vi_end;
  for(tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi)
  {
    std::cout << *vi << std::endl;
  }

  std::cout << "edges: " << std::endl;
  typename graph_traits<Graph_t>::edge_iterator ei, ei_end;
  for(tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
  {
    std::cout << *ei << std::endl;
  }
}

void cndp_cycle_spr::print_spr_vector(vector <set <Vertex> > & my_vector){
	cout<<"============================="<<endl;
	vector <dual_tree_node>::iterator vecIt;

	for (int i = 0; i < my_vector.size(); i++){
			cout<<"................"<<endl;
			cout<<"cycle: ";
			for (set<Vertex>::iterator it2 = my_vector.at(i).begin(); it2 != my_vector.at(i).end(); it2++)
				cout<<*it2<<" ";
			cout<<endl;
	}
	cout<<"******************************"<<endl;
}// end print_dual_tree

void cndp_cycle_spr::print_spr_list(list <fsc_info> & my_list){
	cout<<"============================="<<endl;
	vector <dual_tree_node>::iterator vecIt;
	int a, b, c;
	for (list<fsc_info>::iterator it = my_list.begin(); it != my_list.end(); it++){
		fsc_info cycle = *it;
		a = cycle.out;
		b = cycle.in;
		c = cycle.spr;
//		if ((a >= dim / 3) && (b >= dim / 3)){
			cout<<"................"<<endl;
			cout<<"entrance edge: "<<cycle.entrance_edge<<endl;
			cout<<"cycle: ";
			for (list<Vertex>::iterator it2 = cycle.cycle.begin(); it2 != cycle.cycle.end(); it2++)
				cout<<*it2<<" ";
			cout<<endl;
			cout<<"A: "<<a<<"; B: "<<b<<"; C: "<<c<<" "<<endl;
//  		}
	}
	cout<<"******************************"<<endl;
}// end print_dual_tree

void cndp_cycle_spr::print_dual_tree(vector<dual_tree_node> tree){
	cout<<"+++++++++++++++++++++++++++++++++"<<endl;
	vector <dual_tree_node>::iterator vecIt;
	int a, b, c;
	for (int i = 0; i < tree.size(); i++){
		dual_tree_node dnode = tree.at(i);
		a = dnode.spr_info.out;
		b = dnode.spr_info.in;
		c = dnode.spr_info.spr;

			cout<<"................"<<endl;
			cout<<"Vertex ID: "<<tree.at(i).vertex_id<<endl;
			cout<<"BFS level: "<<dnode.bfs_level<<endl;
			cout<<"Vertex list: "<<tree.at(i).face_vertex_list.at(0)<<" "<<tree.at(i).face_vertex_list.at(1)<<" "<<tree.at(i).face_vertex_list.at(2)<<endl;
			cout<<"edge(0): "<<tree.at(i).face_edge_list.at(0).edge.first<<"-> "<<tree.at(i).face_edge_list.at(0).edge.second<<endl;
			cout<<"edge(1): "<<tree.at(i).face_edge_list.at(1).edge.first<<"-> "<<tree.at(i).face_edge_list.at(1).edge.second<<endl;
			cout<<"edge(2): "<<tree.at(i).face_edge_list.at(2).edge.first<<"-> "<<tree.at(i).face_edge_list.at(2).edge.second<<endl;

			cout<<"entrance edge: "<<dnode.entrance_edge<<endl;

			vector <Edge>::iterator edgeIt;
			cout<<"Outgoing edges: ";
			for (edgeIt = tree.at(i).outgoing_edges.begin(); edgeIt != tree.at(i).outgoing_edges.end(); edgeIt++){
				cout<<*edgeIt<<" ";
			}
			cout<<endl;
			cout<<"child id: "<<" ";
			vector <Vertex>::iterator vtx_it;
			for (vtx_it = tree.at(i).child_id.begin(); vtx_it != tree.at(i).child_id.end(); vtx_it++){
				cout<<*vtx_it<<" ";
			}
			cout<<endl;
			cout<<"cycle: ";
			for (set<Vertex>::iterator it = dnode.f_cycle.begin(); it != dnode.f_cycle.end(); it++)
				cout<<*it<<" ";
			cout<<endl;
			cout<<"A: "<<dnode.spr_info.out<<"; B: "<<dnode.spr_info.in<<"; C: "<<dnode.spr_info.spr<<" "<<endl;

	}
	cout<<"******************************"<<endl;
}// end print_dual_tree

Graph_t cndp_cycle_spr::triangulate(Graph_t &g){
	Graph_t g3 = g;
//	clock_t t = clock();
	//Test for planarity; compute the planar embedding as a side-effect
    typedef std::vector< graph_traits<Graph_t>::edge_descriptor > vec_t;
    std::vector<vec_t> embedding(num_vertices(g3));
    if (boyer_myrvold_planarity_test(boyer_myrvold_params::graph = g3,
                                     boyer_myrvold_params::embedding =
                                         &embedding[0]
                                     )
        )
      std::cout<<"";// "." << std::endl;
    else{
      std::cout << "Input graph is not planar" << std::endl;
      exit(0);
    }

//	clock_t t1 = clock() - t;
//	cout <<"--after planarity test, running time: "<<((float)t1)/CLOCKS_PER_SEC<<endl;


    make_biconnected_planar(g3, &embedding[0]);

//	clock_t t2 = clock() - t;
//	cout <<"--after making biconnected, running time: "<<((float)t2)/CLOCKS_PER_SEC<<endl;

    // Re-initialize the edge index, since we just added a few edges
    graph_traits<Graph_t>::edge_iterator ei, ei_end;
    graph_traits<Graph_t>::edges_size_type edge_count = 0;
    property_map<Graph_t, edge_index_t>::type e_index = get(edge_index, g3);
    edge_count = 0;
    for(boost::tie(ei, ei_end) = edges(g3); ei != ei_end; ++ei)
      put(e_index, *ei, edge_count++);


    //Test for planarity again; compute the planar embedding as a side-effect
    if (boyer_myrvold_planarity_test(boyer_myrvold_params::graph = g3,
                                     boyer_myrvold_params::embedding =
                                         &embedding[0]
                                     )
        )
      std::cout << "";
    else{
      std::cout << "After calling make_biconnected, the graph is not planar"
                << std::endl;
      exit(0);
    }

//	clock_t t3 = clock() - t;
//	cout <<"--after second planarity test, running time: "<<((float)t3)/CLOCKS_PER_SEC<<endl;

    make_maximal_planar(g3, &embedding[0]);
//    clock_t t4 = clock() - t;
//    cout <<"--after making maximal planar, running time: "<<((float)t4)/CLOCKS_PER_SEC<<endl;

return g3;
}
