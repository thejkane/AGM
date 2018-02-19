// Copyright (C) 2004-2012 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Douglas Gregor
//           Nicholas Edmonds
//           Andrew Lumsdaine

#ifndef BOOST_GRAPH_PARALLEL_GRAPHVIZ_HPP
#define BOOST_GRAPH_PARALLEL_GRAPHVIZ_HPP

#ifndef BOOST_GRAPH_USE_MPI
#error "Parallel BGL files should not be included unless <boost/graph/use_mpi.hpp> has been included"
#endif

#include "am++/am++.hpp"

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/distributed/concepts.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/type_traits/is_base_and_derived.hpp>
#include <boost/type_traits/is_same.hpp>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <boost/graph/parallel/container_traits.hpp>
#include <boost/property_map/parallel/global_index_map.hpp>

#include <boost/array.hpp>

struct empty_deleter {void operator()() const {}};

template<typename VertexIndex>
struct msg_vertex_property_data
{
  explicit msg_vertex_property_data() {}
  msg_vertex_property_data(VertexIndex idx, const std::string& p) 
    : idx(idx)
  { strncpy(props, p.c_str(), 64); }

  VertexIndex idx;
  char props[64];
};

struct msg_edge_property_data
{
  explicit msg_edge_property_data() {}
  msg_edge_property_data(const std::string& p)
  { strncpy(props, p.c_str(), 64); }

  char props[64];
};

namespace amplusplus {

  // msg_vertex_property_data
  template <typename A>
  struct make_mpi_datatype<msg_vertex_property_data<A> > : make_mpi_datatype_base {
    make_mpi_datatype<A> dt1;
    scoped_mpi_datatype dt2;
    scoped_mpi_datatype dt;
    make_mpi_datatype(): dt1() {
      int blocklengths[2] = {1, 1};
      MPI_Aint displacements[2];
      msg_vertex_property_data<A> test_object;
      MPI_Aint test_object_ptr;
      MPI_Get_address(&test_object, &test_object_ptr);
      MPI_Get_address(&test_object.idx, &displacements[0]);
      MPI_Get_address(&test_object.props, &displacements[1]);
      displacements[0] -= test_object_ptr;
      displacements[1] -= test_object_ptr;
      MPI_Type_contiguous(64, MPI_CHAR, dt2.get_ptr());
      MPI_Type_commit(dt2.get_ptr());
      MPI_Datatype types[2] = {dt1.get(), dt2.get()};
      MPI_Type_create_struct(2, blocklengths, displacements, types, dt.get_ptr());
      MPI_Type_commit(dt.get_ptr());
    }
    MPI_Datatype get() const {return dt;}
  };

  // msg_edge_property_data
  template <>
  struct make_mpi_datatype<msg_edge_property_data> : make_mpi_datatype_base {
    scoped_mpi_datatype dt1;
    scoped_mpi_datatype dt;
    make_mpi_datatype() {
      int blocklengths[] = {1};
      MPI_Aint displacements[1];
      msg_edge_property_data test_object;
      MPI_Aint test_object_ptr;
      MPI_Get_address(&test_object, &test_object_ptr);
      MPI_Get_address(&test_object.props, &displacements[0]);
      displacements[0] -= test_object_ptr;
      MPI_Type_contiguous(64, MPI_CHAR, dt1.get_ptr());
      MPI_Type_commit(dt1.get_ptr());
      MPI_Datatype types[1] = {dt1.get()};
      MPI_Type_create_struct(1, blocklengths, displacements, types, dt.get_ptr());
      MPI_Type_commit(dt.get_ptr());
    }
    MPI_Datatype get() const {return dt;}
  };
} // namespace amplusplus

template<typename VerticesSizeType>
struct local_sizes_handler {

  local_sizes_handler() {}
  local_sizes_handler(std::vector<VerticesSizeType>& sizes) : sizes(&sizes) {}

  void operator() (amplusplus::transport::rank_type source, 
		   const VerticesSizeType* local_size, size_t /* count */) const
  {
    (*sizes)[source] = *local_size;
  }

protected:
  std::vector<VerticesSizeType>* sizes;
};

template<typename ProcessSizeType>
struct vertex_property_handler {

  vertex_property_handler() {}
  vertex_property_handler(std::vector<std::pair<ProcessSizeType, std::string> >& vp) 
    : vertex_props(&vp) {}

  template<typename VertexIndex>
  void operator() (int source, const msg_vertex_property_data<VertexIndex>& data)
  {
    assert(data.idx < vertex_props->size());
    (*vertex_props)[data.idx] = std::make_pair(source, data.props);
  }

protected:
  std::vector<std::pair<ProcessSizeType, std::string> >* vertex_props;
};

struct edge_property_handler {

  edge_property_handler() {}
  edge_property_handler(std::vector<std::vector<std::string> >& ep) 
    : edge_props(&ep) {}

  void operator() (int source, const msg_edge_property_data& data)
  {
    (*edge_props)[source].push_back(data.props);
  }

protected:
  std::vector<std::vector<std::string> >* edge_props;
};

namespace boost {

template<typename Graph>
struct graph_id_writer
{
  explicit graph_id_writer(Graph& g) : g(g) { }

  void operator()(std::ostream& out)
  {
    out << "    label=\"p" << g->process_group().rank() << "\";\n";
  }

 private:
  Graph& g;
};

template<typename NumberMap>
struct paint_by_number_writer
{
  explicit paint_by_number_writer(NumberMap number) : number(number) { }

  template<typename Descriptor>
  void operator()(std::ostream& out, Descriptor k)
  {
    static const char* color_names[] = {
      "blue",
      "brown",
      "cyan",
      "darkgreen",
      "darkorchid",
      "darksalmon",
      "darkviolet",
      "deeppink",
      "gold3",
      "green",
      "magenta",
      "navy",
      "red",
      "yellow",
      "palegreen",
      "gray65",
      "gray21",
      "bisque2",
      "greenyellow",
      "indianred4",
      "lightblue2",
      "mediumspringgreen",
      "orangered",
      "orange"
    };
    const int colors = sizeof(color_names) / sizeof(color_names[0]);
    if (get(number, k) < colors) {
      out << " [ style=\"filled\", fillcolor=\"" << color_names[get(number, k)]
          << "\" ]";
    } else {
      out << " [ label=\"(" << get(number, k) << ")\" ]";
    }
  }

 private:
  NumberMap number;
};

template<typename NumberMap>
inline paint_by_number_writer<NumberMap>
paint_by_number(NumberMap number)
{ return paint_by_number_writer<NumberMap>(number); }

template<typename Graph, typename VertexPropertiesWriter, 
         typename EdgePropertiesWriter, typename GraphPropertiesWriter>
void 
write_graphviz(std::ostream& out,
               Graph& g, 
               VertexPropertiesWriter vpw,
               EdgePropertiesWriter epw,
               GraphPropertiesWriter /* gpw unused */
               BOOST_GRAPH_ENABLE_IF_MODELS_PARM(Graph,distributed_graph_tag))
{
  typedef typename graph_traits<Graph>::directed_category directed_category;
  typedef typename graph_traits<Graph>::vertices_size_type vertices_size_type;
  typedef typename graph_traits<Graph>::edges_size_type edges_size_type;
  typedef typename amplusplus::transport::rank_type rank_type;
  typedef typename amplusplus::transport::rank_type process_size_type;
  typedef typename property_map<Graph, vertex_index_t>::type VertexIndexMap;
  typedef typename property_map<Graph, vertex_global_t>::type VertexGlobalMap;

  rank_type id = g.transport().rank(); 
  process_size_type N = g.transport().size(); 

  static const bool is_undirected
    = (is_base_and_derived<undirected_tag, directed_category>::value
       || is_same<undirected_tag, directed_category>::value);
  static const char* graph_kind = is_undirected? "graph" : "digraph";
  static const char* edge_kind = is_undirected? "--" : "->";

  parallel::global_index_map<VertexIndexMap, VertexGlobalMap> 
    global_index(g.transport(), num_vertices(g), get(vertex_index, g),
                 get(vertex_global, g));
  // Find total number of vertices
  // TODO: This would be much easier with all_reduce()
  std::vector<vertices_size_type> n_verts(N, 0);
  std::vector<edges_size_type> n_edges(N, 0);
  n_verts[id] = num_vertices(g);
  n_edges[id] = num_edges(g);

  // Setup messages for sizes... can we use all_reduce here?
  amplusplus::message_type<vertices_size_type> vertices_size_msg = 
    g.transport().template create_message_type<vertices_size_type>();
  amplusplus::message_type<edges_size_type> edges_size_msg = 
    g.transport().template create_message_type<edges_size_type>();

  vertices_size_msg.set_max_count(1);
  edges_size_msg.set_max_count(1);

  vertices_size_msg.set_handler(local_sizes_handler<vertices_size_type>(n_verts));
  edges_size_msg.set_handler(local_sizes_handler<edges_size_type>(n_edges));
  {
    amplusplus::scoped_epoch epoch(g.transport());

    vertices_size_type v = num_vertices(g);
    edges_size_type e = num_edges(g);
    if (id != 0) {
      vertices_size_msg.message_being_built(0);
      edges_size_msg.message_being_built(0);
      vertices_size_msg.send(&v, 1, 0, empty_deleter());
      edges_size_msg.send(&e, 1, 0, empty_deleter());
    }
  }

  // Allocate space for vertex properties and reserve space for edge properties
  std::vector<std::pair<rank_type, std::string> > vertex_property_data(N);
  std::vector<std::vector<std::string> > edge_property_data(N);
  if (id == 0) {
    vertex_property_data.resize(std::accumulate(n_verts.begin(), n_verts.end(), 0, 
						std::plus<vertices_size_type>()));
    for (process_size_type i = 0; i < N; ++i)
      edge_property_data[i].reserve(n_edges[i]);
  }

  // Setup vertex property messages
  typedef msg_vertex_property_data<vertices_size_type> msg_vp_data;

  typedef amplusplus::basic_coalesced_message_type<msg_vp_data,
                                                   vertex_property_handler<rank_type> >
    vertex_property_message_type;

  amplusplus::register_mpi_datatype<msg_vp_data>();
  vertex_property_message_type vertex_property_msg(amplusplus::basic_coalesced_message_type_gen(1 << 12), g.transport());
  vertex_property_msg.set_handler(vertex_property_handler<rank_type>(vertex_property_data));

  // Setup edge property messages
  typedef amplusplus::basic_coalesced_message_type<msg_edge_property_data,
                                                   edge_property_handler>
    edge_property_message_type;

  amplusplus::register_mpi_datatype<msg_edge_property_data>();
  edge_property_message_type edge_property_msg(amplusplus::basic_coalesced_message_type_gen(1 << 12), g.transport());
  edge_property_msg.set_handler(edge_property_handler(edge_property_data));

  // Send properties to process 0
  {
    amplusplus::scoped_epoch epoch(g.transport());

    // TODO (NGE): Send and write graph properties
    
    // Vertex properties
    typename graph_traits<Graph>::vertex_iterator vi, vi_end;
    for (tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi) {
      std::ostringstream vertex_prop;
      vpw(vertex_prop, *vi);
      vertex_property_msg.send(msg_vp_data(get(global_index, *vi), vertex_prop.str()), 0);
    }
    
    // Edge properties
    typename graph_traits<Graph>::edge_iterator ei, ei_end;
    for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
      std::ostringstream edge_prop;
      int source_idx = get(global_index, source(*ei, g));
      int target_idx = get(global_index, target(*ei, g));
      edge_prop << "  n" << source_idx << " " << edge_kind << " n" << target_idx;
      epw(edge_prop, *ei);
      edge_prop << ";\n";
      edge_property_msg.send(msg_edge_property_data(edge_prop.str()), 0);
    }
  }

  // Receive properties on rank 0 and write output file
  if (id == 0) {
    out << graph_kind << " g {\n";
    for (rank_type id = 0; id < N; ++id) {
      out << "  subgraph cluster_" << id << " {\n";
      for (std::size_t i = 0; i < vertex_property_data.size(); ++i)
	if (vertex_property_data[i].first == id)
	  out << "    n" << i << vertex_property_data[i].second << ";\n";
      out << "  }\n\n";

      // Write out edge properties, edges should be ordered by channel semantics
      for (edges_size_type i = 0; i < edge_property_data[id].size(); ++i)
	out << edge_property_data[id][i];
    }
    out << "}\n";
  }
}

template<typename Graph, typename VertexPropertiesWriter, 
         typename EdgePropertiesWriter>
inline void 
write_graphviz(std::ostream& out,
               Graph& g, 
               VertexPropertiesWriter vpw,
               EdgePropertiesWriter epw
               BOOST_GRAPH_ENABLE_IF_MODELS_PARM(Graph,distributed_graph_tag))
{
  write_graphviz(out, g, vpw, epw, graph_id_writer<Graph>(g));
}

template<typename Graph, typename VertexPropertiesWriter>
inline void 
write_graphviz(std::ostream& out,
               Graph& g, 
               VertexPropertiesWriter vpw
               BOOST_GRAPH_ENABLE_IF_MODELS_PARM(Graph,distributed_graph_tag))
{
  write_graphviz(out, g, vpw, default_writer());
}

template<typename Graph>
inline void 
write_graphviz(std::ostream& out, Graph& g
               BOOST_GRAPH_ENABLE_IF_MODELS_PARM(Graph,distributed_graph_tag))
{
  write_graphviz(out, g, default_writer());
}

template<typename Graph, typename VertexPropertiesWriter, 
         typename EdgePropertiesWriter, typename GraphPropertiesWriter>
void 
write_graphviz(const std::string& filename,
               Graph& g, 
               VertexPropertiesWriter vpw,
               EdgePropertiesWriter epw,
               GraphPropertiesWriter gpw
               BOOST_GRAPH_ENABLE_IF_MODELS_PARM(Graph,distributed_graph_tag))
{
  if (g.transport().rank() == 0) {
    std::ofstream out(filename.c_str());
    write_graphviz(out, g, vpw, epw, gpw);
  } else {
    write_graphviz(std::cout, g, vpw, epw, gpw);
  }
}

template<typename Graph, typename VertexPropertiesWriter, 
         typename EdgePropertiesWriter>
void 
write_graphviz(const std::string& filename,
               Graph& g, 
               VertexPropertiesWriter vpw,
               EdgePropertiesWriter epw
               BOOST_GRAPH_ENABLE_IF_MODELS_PARM(Graph,distributed_graph_tag))
{
  if (g.transport().rank() == 0) {
    std::ofstream out(filename.c_str());
    write_graphviz(out, g, vpw, epw);
  } else {
    write_graphviz(std::cout, g, vpw, epw);
  }
}

template<typename Graph, typename VertexPropertiesWriter>
void 
write_graphviz(const std::string& filename,
               Graph& g, 
               VertexPropertiesWriter vpw
               BOOST_GRAPH_ENABLE_IF_MODELS_PARM(Graph,distributed_graph_tag))
{
  if (g.transport().rank() == 0) {
    std::ofstream out(filename.c_str());
    write_graphviz(out, g, vpw);
  } else {
    write_graphviz(std::cout, g, vpw);
  }
}

template<typename Graph>
void 
write_graphviz(const std::string& filename, Graph& g
               BOOST_GRAPH_ENABLE_IF_MODELS_PARM(Graph,distributed_graph_tag))
{
  if (g.transport().rank() == 0) {
    std::ofstream out(filename.c_str());
    write_graphviz(out, g);
  } else {
    write_graphviz(std::cout, g);
  }
}

template<typename Graph>
void
write_graphviz(std::ostream& out, Graph& g,
               const dynamic_properties& dp, 
               const std::string& node_id = "node_id"
               BOOST_GRAPH_ENABLE_IF_MODELS_PARM(Graph,distributed_graph_tag))
{
  write_graphviz
    (out, g,
     /*vertex_writer=*/dynamic_vertex_properties_writer(dp, node_id),
     /*edge_writer=*/dynamic_properties_writer(dp));
}

} // end namespace boost

#endif // BOOST_GRAPH_PARALLEL_GRAPHVIZ_HPP
