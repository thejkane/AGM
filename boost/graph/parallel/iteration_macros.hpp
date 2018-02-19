// Copyright (C) 2011-2012 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Nicholas Edmonds
//           Thejaka Kanewala
//           Andrew Lumsdaine

#ifndef BOOST_GRAPH_PARALLEL_ITERATION_MACROS_HPP
#define BOOST_GRAPH_PARALLEL_ITERATION_MACROS_HPP

#include <utility>
#include <boost/graph/iteration_macros.hpp>

#define BGL_VSTART(linenum) BGL_CAT(bgl_vstart_,linenum)
#define BGL_VEND(linenum) BGL_CAT(bgl_vend_,linenum)
#define BGL_ESTART(linenum) BGL_CAT(bgl_estart_,linenum)
#define BGL_EEND(linenum) BGL_CAT(bgl_eend_,linenum)

#define BGL_PARFORALL_VERTICES(VNAME, GNAME, GraphType, TID, NTHREADS) \
for (boost::graph_traits<GraphType>::vertex_iterator \
     BGL_VSTART(__LINE__) \
       = (vertices(GNAME).first + (TID * (num_vertices(GNAME) + NTHREADS) / NTHREADS)), \
     BGL_VEND(__LINE__) \
       = (vertices(GNAME).first + std::min BOOST_PREVENT_MACRO_SUBSTITUTION \
           ((TID + 1) * (num_vertices(GNAME) + NTHREADS) / NTHREADS, num_vertices(GNAME))); \
  BGL_VSTART(__LINE__) != BGL_VEND(__LINE__); BGL_VSTART(__LINE__) = BGL_VEND(__LINE__)) \
  for (boost::graph_traits<GraphType>::vertex_descriptor VNAME; \
    BGL_VSTART(__LINE__) != BGL_VEND(__LINE__) ? (VNAME = *BGL_VSTART(__LINE__), true):false; \
     ++BGL_VSTART(__LINE__))

// This only works for RandomAccess Iterators
#define BGL_PARFORALL_VERTICES_T(VNAME, GNAME, GraphType, TID, NTHREADS) \
for (typename boost::graph_traits<GraphType>::vertex_iterator \
     BGL_VSTART(__LINE__) \
       = (vertices(GNAME).first + (TID * (num_vertices(GNAME) + NTHREADS) / NTHREADS)), \
     BGL_VEND(__LINE__) \
       = (std::min)(vertices(GNAME).first + ((TID + 1) * (num_vertices(GNAME) + NTHREADS) / NTHREADS), vertices(GNAME).second); \
  BGL_VSTART(__LINE__) != BGL_VEND(__LINE__); BGL_VSTART(__LINE__) = BGL_VEND(__LINE__)) \
  for (typename boost::graph_traits<GraphType>::vertex_descriptor VNAME; \
    BGL_VSTART(__LINE__) < BGL_VEND(__LINE__) ? (VNAME = *BGL_VSTART(__LINE__), true):false; \
     ++BGL_VSTART(__LINE__))

// This only works for RandomAccess Iterators
#define BGL_PARFORALL_EDGES_T(ENAME, GNAME, GraphType, TID, NTHREADS) \
for (typename boost::graph_traits<GraphType>::edge_iterator \
     BGL_ESTART(__LINE__) \
       = (edges(GNAME).first + (TID * (num_edges(GNAME) + NTHREADS) / NTHREADS)), \
     BGL_EEND(__LINE__) \
       = (std::min)(edges(GNAME).first + ((TID + 1) * (num_edges(GNAME) + NTHREADS) / NTHREADS), edges(GNAME).second); \
  BGL_ESTART(__LINE__) != BGL_EEND(__LINE__); BGL_ESTART(__LINE__) = BGL_EEND(__LINE__)) \
  for (typename boost::graph_traits<GraphType>::edge_descriptor ENAME; \
    BGL_ESTART(__LINE__) < BGL_EEND(__LINE__) ? (ENAME = *BGL_ESTART(__LINE__), true):false; \
     ++BGL_ESTART(__LINE__))

#endif // BOOST_GRAPH_PARALLEL_ITERATION_MACROS_HPP
