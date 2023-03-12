#ifndef __CREATE_DUAL_GRAPH_HPP__
#define __CREATE_DUAL_GRAPH_HPP__
#include <tuple>
#include <vector>
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/planar_face_traversal.hpp>

namespace boost
{
	
  template <typename InputGraph,
            typename OutputGraph,
            typename EdgeIndexMap,
			typename EdgeWeightMap
			>
  struct dual_graph_visitor : public planar_face_traversal_visitor
  {
    typedef typename graph_traits<OutputGraph>::vertex_descriptor vertex_t;
    typedef typename graph_traits<InputGraph>::edge_descriptor edge_t;
    typedef typename std::vector<vertex_t> vertex_vector_t;
    typedef iterator_property_map< typename vertex_vector_t::iterator, EdgeIndexMap >    edge_to_face_map_t;
    typedef typename std::vector <std::vector <vertex_t> > vertex_vector_vector_t;



//    typedef typename std::vector<int> vector_t;
	// declaring properties
    InputGraph& g;
    OutputGraph& dual_g;
    EdgeIndexMap em;
	EdgeWeightMap wm;
    vertex_t current_face;
    vertex_vector_t edge_to_face_vector;
    edge_to_face_map_t edge_to_face;
	vertex_vector_vector_t faceVector;
	vertex_vector_t tmpVector;

	// constructor	
    dual_graph_visitor(InputGraph& arg_g,
                       OutputGraph& arg_dual_g,
                       EdgeIndexMap arg_em,
					   EdgeWeightMap arg_wm
                       ) :
      g(arg_g),
      dual_g(arg_dual_g),
      em(arg_em),
	  wm(arg_wm),
      edge_to_face_vector(num_edges(g),
                          graph_traits<OutputGraph>::null_vertex()),
      edge_to_face(edge_to_face_vector.begin(), em)
      {}


	
  
// re-defining the event in planar_face_traversal_visitor class
    void begin_traversal()
    {
//		std::cout<<"~~~Hello from planar_dual_tree.begin_traversal(): "<<std::endl;				
	}
	
// re-defining the event in planar_face_traversal_visitor class	
    void begin_face()
    {
//		std::cout<<"New face: ";
		current_face = add_vertex(dual_g);
		tmpVector.resize(0);// = vertex_vector_t;
		//		vertex_counter = 0;
    }


// re-defining the event in planar_face_traversal_visitor class
   template <typename Vertex>
   void next_vertex(Vertex v)
   {
//   		std::cout << v << " ";
   		tmpVector.push_back(v);
//   		vertex_counter++;
   }

// re-defining the event in planar_face_traversal_visitor class
	void end_face(){
//		std::cout << std::endl;
		faceVector.push_back(tmpVector);
	}

// re-defining the event in planar_face_traversal_visitor class
    template <typename Edge>
    void next_edge(Edge e)
    {
        typedef typename property_map <InputGraph, edge_weight_t>::type pm_edgeWeightType;
    	typedef typename boost::property_traits<pm_edgeWeightType>::value_type pt_edgeWeightType;
    	pt_edgeWeightType my_weight;
      vertex_t existing_face = edge_to_face[e];
      if (existing_face == graph_traits<OutputGraph>::null_vertex())
      {
        edge_to_face[e] = current_face;
      }
      else
      {
    	  my_weight = boost::get(edge_weight, g, e);
    	  if (my_weight == 0){
    		  add_edge(existing_face, current_face, dual_g);
    	  }
      }
    }// end next_edge()

// defining getting face_vector
	vertex_vector_vector_t get_face_vector(){
		return faceVector;
	}
		
};// end dual_graph_visitor class


  
// create_dual_graph function with 4 arguments
  template <typename InputGraph,
            typename OutputGraph,
            typename PlanarEmbedding,
            typename EdgeIndexMap,
			typename EdgeWeightMap				
           >
  std::vector < std::vector <typename graph_traits<OutputGraph>::vertex_descriptor> > 
  create_dual_graph(InputGraph& g,
                         OutputGraph& dual_g,
                         PlanarEmbedding embedding,
                         EdgeIndexMap em,
						 EdgeWeightMap wm
						 )
  {
    dual_graph_visitor<InputGraph, OutputGraph, EdgeIndexMap, EdgeWeightMap>
	visitor(g, dual_g, em, wm);
    planar_face_traversal(g, embedding, visitor, em);
	typedef typename graph_traits<OutputGraph>::vertex_descriptor vertex_t;
	typedef typename graph_traits<OutputGraph>::vertex_iterator vertex_iter_t;
    std::vector <std::vector <vertex_t> > faceVector;
    typedef typename std::vector <std::vector <vertex_t> >::iterator fv_iterator_t;
	typedef typename std::vector <vertex_t>::iterator vertex_iterator_t;

	fv_iterator_t fv_iter;
	vertex_iterator_t vi;
	faceVector = visitor.get_face_vector();
	// print out the faceVector
//	for (fv_iter = faceVector.begin(); fv_iter != faceVector.end(); fv_iter++){
//		for (vi = fv_iter->begin(); vi != fv_iter->end(); vi++){
//			std::cout<<*vi<<" ";
//		}
//		std::cout<<std::endl;
//	}

	return faceVector;
  }// end create_dual_graph()

// create_dual_graph function with 3 arguments
/* 
 template <typename InputGraph,
            typename OutputGraph,
            typename PlanarEmbedding
           >
  void create_dual_graph(InputGraph& g,
                         OutputGraph& dual_g,
                         PlanarEmbedding embedding
                         )
  {
    create_dual_graph(g, dual_g, embedding, get(edge_index,g), get(edge_weight,g));
  }*/
} // namespace boost

#endif //__CREATE_DUAL_GRAPH_HPP__
