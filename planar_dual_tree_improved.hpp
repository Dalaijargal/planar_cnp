#ifndef __CREATE_DUAL_GRAPH_HPP__
#define __CREATE_DUAL_GRAPH_HPP__

#include <vector>

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/planar_face_traversal_improved.hpp>

namespace boost
{

  template <typename InputGraph,
            typename OutputGraph,
            typename EdgeIndexMap,
			typename EdgeWeightMap,
			typename Vertex
			>
  struct dual_graph_visitor : public planar_face_traversal_visitor
  {

    typedef typename graph_traits<OutputGraph>::vertex_descriptor vertex_t;
    typedef typename graph_traits<InputGraph>::edge_descriptor edge_t;
    typedef typename std::vector<vertex_t> vertex_vector_t;
    typedef iterator_property_map
      < typename vertex_vector_t::iterator, EdgeIndexMap >
        edge_to_face_map_t;

    dual_graph_visitor(InputGraph& arg_g,
                       OutputGraph& arg_dual_g,
                       EdgeIndexMap arg_em,
					   EdgeWeightMap arg_wm,
					   Vertex arg_v1,
					   Vertex arg_v2,
					   Vertex arg_v3
                       ) :
      g(arg_g),
      dual_g(arg_dual_g),
      em(arg_em),
	  wm(arg_wm),
	  v1(arg_v1),
	  v2(arg_v2),
	  v3(arg_v3),
      edge_to_face_vector(num_edges(g),
                          graph_traits<OutputGraph>::null_vertex()),
      edge_to_face(edge_to_face_vector.begin(), em)
      {}

//	template <typename Vertex>
    void begin_traversal(Vertex v1, Vertex v2, Vertex v3)
    {
		std::cout<<"Hello from begin_traversal: "<<v1<<" "<<v2<<" "<<v3<<std::endl;				
	}
	
    void begin_face()
    {
		current_face = add_vertex(dual_g);
    }

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
//		  std::cout<<"~***"<<my_weight<<std::endl;
    	  if (my_weight == 0){
    		  add_edge(existing_face, current_face, dual_g);
    		  //std::cout<<"~"<<my_weight<<std::endl;
    	  }
      }
    }

    InputGraph& g;
    OutputGraph& dual_g;
    EdgeIndexMap em;
	EdgeWeightMap wm;
	Vertex v1;
	Vertex v2;
	Vertex v3;
    vertex_t current_face;
    vertex_vector_t edge_to_face_vector;
    edge_to_face_map_t edge_to_face;
  };

  template <typename InputGraph,
            typename OutputGraph,
            typename PlanarEmbedding,
            typename EdgeIndexMap,
			typename EdgeWeightMap,
			typename Vertex			
           >
  void create_dual_graph(InputGraph& g,
                         OutputGraph& dual_g,
                         PlanarEmbedding embedding,
                         EdgeIndexMap em,
						 EdgeWeightMap wm,
						 Vertex v1,
						 Vertex v2,
						 Vertex v3)
  {
    dual_graph_visitor<InputGraph, OutputGraph, EdgeIndexMap, EdgeWeightMap, Vertex>   
	visitor(g, dual_g, em, wm, v1, v2, v3);
					   
    planar_face_traversal(g, embedding, visitor, em, v1, v2, v3);
  }

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
  }

} // namespace boost

#endif //__CREATE_DUAL_GRAPH_HPP__
