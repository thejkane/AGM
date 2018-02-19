// Copyright (C) 2007-2012 The Trustees of Indiana University.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Thejaka Kanewala
//           Marcin Zalewski
//           Andrew Lumsdaine
//======== Connected Components Algortihm================//
// Delta-Stepping version of connected components algorithm.
//===========================================================//


#ifndef BOOST_GRAPH_DELTA_STEPPING_CC_HPP
#define BOOST_GRAPH_DELTA_STEPPING_CC_HPP

#ifndef BOOST_GRAPH_USE_MPI
#error "Parallel BGL files should not be included unless <boost/graph/use_mpi.hpp> has been included"
#endif

#include <am++/counter_coalesced_message_type.hpp>
#include <am++/detail/thread_support.hpp>

#include <boost/parallel/append_buffer.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/parallel/algorithm.hpp> // for all_reduce
#include <boost/graph/parallel/thread_support.hpp> // for compare_and_swap
#include <algorithm> // for std::min, std::max
#include <boost/format.hpp>
#include <iostream>
#include <atomic>

#include <map>


namespace boost { namespace graph { namespace distributed {

template<typename Graph, 
	 typename ComponentMap, 
	 typename IdDistribution,
	 typename WorkStats,
         typename MessageGenerator = 
           amplusplus::simple_generator<amplusplus::counter_coalesced_message_type_gen> >
class delta_stepping_cc {

  typedef delta_stepping_cc<Graph, ComponentMap, IdDistribution, WorkStats, MessageGenerator> 
  self_type;

  typedef typename boost::property_map<Graph, vertex_owner_t>::const_type OwnerMap;
  typedef typename property_traits<ComponentMap>::value_type Component;
  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  typedef typename graph_traits<Graph>::degree_size_type Degree;
  typedef std::pair<Vertex, Component> vertex_component_data;

  typedef append_buffer<vertex_component_data, 10u> Bucket;
  typedef append_buffer<vertex_component_data, 10u> InitialBuffer;

  typedef typename graph_traits<Graph>::vertices_size_type VerticesSize;
  typedef typename std::vector<Bucket*>::size_type BucketIndex;
  typedef typename graph_traits<Graph>::degree_size_type DegreeSize;
  
  struct vertex_component_handler;

  typedef typename MessageGenerator::template call_result<vertex_component_data, vertex_component_handler, 
							  vertex_component_owner<OwnerMap, vertex_component_data>,
							  amplusplus::idempotent_combination_t<boost::parallel::minimum<Component>,
											       Component> >::type
  RelaxMessage;

public:
  delta_stepping_cc(Graph& g,
		    amplusplus::transport &t,
		    int offs,
		    ComponentMap components, 
		    BucketIndex nbkts,
		    IdDistribution& idd,
		    WorkStats& stats,
		    MessageGenerator message_gen = 
		    MessageGenerator(amplusplus::counter_coalesced_message_type_gen(1 << 12)))
    : dummy_first_member_for_init_order((amplusplus::register_mpi_datatype<vertex_component_data>(), 0)),
      g(g), 
      transport(t), 
      core_offset(offs),
      components(components), 
      num_buckets(nbkts), 
      owner(get(vertex_owner, g)), 
      current_level(0), 
    buckets_processed(0), current_bucket(0), 
    id_distribution(idd),
    work_stats(stats),
    relax_msg(message_gen, transport, vertex_component_owner<OwnerMap, vertex_component_data>(owner),
	      amplusplus::idempotent_combination(boost::parallel::minimum<Component>(), 
						 std::numeric_limits<Component>::max()))
  {
    initialize();
  }

  delta_stepping_cc(Graph& g,
		    amplusplus::transport &t,
		    int offs,
		    ComponentMap components, 
		    IdDistribution& idd,
		    WorkStats& stats,
		    MessageGenerator message_gen = 
		    MessageGenerator(amplusplus::counter_coalesced_message_type_gen(1 << 12)))
    : dummy_first_member_for_init_order((amplusplus::register_mpi_datatype<vertex_component_data>(), 0)),
      g(g), 
      transport(t), 
      core_offset(offs),
      components(components), delta(0), num_buckets(0), owner(get(vertex_owner, g)), 
      current_level(0), buckets_processed(0), current_bucket(0), 
    id_distribution(idd),
    work_stats(stats),
    relax_msg(message_gen, transport, vertex_component_owner<OwnerMap, vertex_component_data>(owner),
	      amplusplus::idempotent_combination(boost::parallel::minimum<Component>(), 
						 std::numeric_limits<Component>::max()))
  {
    initialize();
  }



#ifdef AMPLUSPLUS_PRINT_HIT_RATES
  typedef std::pair<unsigned long long, unsigned long long> cache_stats;
  cache_stats
  get_cache_stats() {
    cache_stats stats(0, 0);
    for(size_t i = 0; i < relax_msg.counters.hits.size(); ++i) {
      stats.first += relax_msg.counters.hits[i];
      stats.second += relax_msg.counters.tests[i];
    }
    return stats;
  }
#endif // AMPLUSPLUS_PRINT_HIT_RATES

  void operator() (int tid) { run(tid); }

  void run(int tid = 0);
  void populate_initial_buffer(int tid);

  BucketIndex get_num_levels() { return buckets_processed; }
  time_type get_start_time() { return start_time; }
  time_type get_end_time() { return end_time; }

  template<typename SizeType>
  SizeType logical_id(SizeType k) {
    return id_distribution(k);
  }


protected:
  void initialize();
  
  // Relax the edge (u, v), creating a new best path of components x.
  void relax(const vertex_component_data& vdc,
	     int tid);

  void find_next_bucket();

  const int dummy_first_member_for_init_order; // Unused

  const Graph& g;
  amplusplus::transport& transport;
  int core_offset;
  ComponentMap components;
  VerticesSize delta;
  OwnerMap owner;

  // Bucket data structure. The ith bucket contains all local vertices
  // with (tentative) components in the range [i*delta,
  // (i+1)*delta). 
  std::vector<shared_ptr<Bucket> > buckets;

  // Bucket to hold vertices deleted at each level
  BucketIndex current_level; // How many buckets have we processed?
  BucketIndex buckets_processed; // Stats tracking
  BucketIndex num_buckets;

  // Shared thread state to make sure we're all on the same page
  BucketIndex current_bucket;
  shared_ptr<amplusplus::detail::barrier> t_bar;
  RelaxMessage relax_msg;
  time_type start_time;
  time_type end_time;
  IdDistribution& id_distribution;
  shared_ptr<InitialBuffer> buffer;
  WorkStats& work_stats;
};

#define DELTA_STEPPING_CC_PARMS                                   \
      typename Graph, typename ComponentMap, typename IdDistribution, typename WorkStats, typename MessageGenerator

#define DELTA_STEPPING_CC_TYPE                                    \
      delta_stepping_cc<Graph, ComponentMap, IdDistribution, WorkStats, MessageGenerator>

template<DELTA_STEPPING_CC_PARMS>
void
DELTA_STEPPING_CC_TYPE::initialize() {
  relax_msg.set_handler(vertex_component_handler(*this));

  // calculate initial set
  // also calculate the max vertex in the initial
  // vertex set
  Vertex localmaxv = 0;

  { 
    BGL_FORALL_VERTICES_T(u, g, Graph) {
      Vertex min_neighbor = logical_id(u);
      std::set<Vertex> adjacencies1;
      BGL_FORALL_OUTEDGES_T(u, e, g, Graph) {
	Vertex lv = logical_id(target(e, g));
	if (logical_id(u) != lv) { // ignore self loops
	  if (adjacencies1.insert(lv).second) {
	    // check whether u is the minimum out of its neighbors
	    if (lv < min_neighbor) {
	      min_neighbor = lv;
	      break;
	    }
	  }
	}
      }

      if (min_neighbor == logical_id(u)) {
	if (logical_id(u) > localmaxv)
	  localmaxv = logical_id(u);
      }
    }
  } // end of epoch 

  //  std::cout << "Local maximum vertex : " << localmaxv << std::endl;

  // find the global maxv

  using boost::parallel::all_reduce;
  using boost::parallel::maximum;
  using std::max;

  all_reduce<Vertex, maximum<Vertex> > r(transport, maximum<Vertex>());
  Vertex maxv = r(localmaxv);

  if (transport.rank() == 0)
    std::cout << "Maximum vertex id in the initial set : " << maxv << std::endl;

  // If num buckets wasn't supplied in the ctor initialize it
  if (num_buckets == 0) {
    // If num buckets is not supplied take a guess at
    // num buckets. This is not optimal just a guess
    num_buckets = 6; 
    // why 6 ? well in most cases graphs has low
    // diameters therefore assume there will be 6 buckets.
    // TODO change later
  }

  delta = maxv / num_buckets;

  if (delta == 0) {
    std::cout << " num_buckets > maxv, num_buckets = " << num_buckets
	      << ", maxv = " << maxv << " resetting buckets to 1" <<  std::endl;
    delta = maxv;
  }

  // need an additional bucket
  num_buckets += 1;
  
  if (transport.rank() == 0)
    std::cout << "Delta value = " << delta << std::endl;

  // Declare bucket data structure and index variable
  buckets.resize(num_buckets);
  for (BucketIndex i = 0 ; i < buckets.size() ; ++i) {
    shared_ptr<Bucket> p(new Bucket);
    buckets[i].swap(p);
  }

  // Initialize components labels
  BGL_FORALL_VERTICES_T(v, g, Graph) { 
    put(components, v, logical_id(v)); 
  }

  // Initial buffer
  shared_ptr<InitialBuffer> p(new InitialBuffer);
  buffer.swap(p);
}


template<DELTA_STEPPING_CC_PARMS>
void
DELTA_STEPPING_CC_TYPE::populate_initial_buffer(int tid) {

  int nthreads = transport.get_nthreads();

  BGL_PARFORALL_VERTICES_T(u, g, Graph, tid, nthreads) {
    Vertex min_neighbor = logical_id(u);
    std::set<Vertex> adjacencies1;
    BGL_FORALL_OUTEDGES_T(u, e, g, Graph) {
      Vertex v = target(e, g);
      if (u != v) { // ignore self loops
	if (adjacencies1.insert(v).second) {
	  // check whether u is the minimum out of its neighbors
	  if (logical_id(v) < min_neighbor) {
	    min_neighbor = logical_id(v);
	    break;
	  }
	}
      }
    }

    if (min_neighbor == logical_id(u)) {
      std::set<Vertex> adjacencies;
      BGL_FORALL_OUTEDGES_T(u, e, g, Graph) {
	Vertex v = target(e, g);
	if (v != u) {
	  if (adjacencies.insert(v).second) {
	    buffer->push_back(vertex_component_data(v, logical_id(u)));
	  }
	}
      }
    }
  }
}

template<DELTA_STEPPING_CC_PARMS>
void
DELTA_STEPPING_CC_TYPE::run(int tid) {

  int count_epoch = 0 ; 
  AMPLUSPLUS_WITH_THREAD_ID(tid) {
    int nthreads = transport.get_nthreads();
    
    if (tid == 0)
      t_bar.reset(new amplusplus::detail::barrier(nthreads));
    
    // This barrier acts as a temporary barrier until we can be sure t_bar is initialized 
    { amplusplus::scoped_epoch epoch(transport); } 

    t_bar->wait();
    if (pin(tid+core_offset) != 0) {
      std::cerr << "[ERROR] Unable to pin current thread to "
		<< "core : " << tid << std::endl;
      assert(false);
    }

    // wait till all threads are pinned
    t_bar->wait();

    validate_thread_core_relation();

    // Now above if branch needs to be executed to every thread
    // Therefore wait till every thread comes to this point
    t_bar->wait();

    // last passing thread will update the correct start time
    populate_initial_buffer(tid);

    { amplusplus::scoped_epoch epoch(transport); }

    typename InitialBuffer::size_type buff_size;
    buff_size = buffer->size();

    // start timer
    start_time = get_time();

    { 
      amplusplus::scoped_epoch epoch(transport);
      for (typename InitialBuffer::size_type i = tid ; 
	   i < buff_size ; i+= nthreads) {
	vertex_component_data& vd = (*buffer)[i];
	relax_msg.send(vd);
      }
    }

    // clear the initial buffer
    buffer->clear();

    t_bar->wait();
    
    typename Bucket::size_type current_bucket_start, current_bucket_end;
    
    do {
      
      current_bucket_start = 0;
      current_bucket_end = buckets[current_bucket]->size();
      
      // Process current bucket
      do {
        
        unsigned long all_starting_sizes;
        const unsigned long starting_size = current_bucket_end - current_bucket_start;
        count_epoch = 0;
        t_bar->wait();
        {
          amplusplus::scoped_epoch_value epoch(transport, starting_size, all_starting_sizes);
          
#ifdef PRINT_DEBUG
          if (tid == 0)
            std::cout << transport.rank() << ": processing bucket "
                      << current_bucket << " [" << current_bucket_start << ", " << current_bucket_end << ")\n";
#endif
	  
	  while (current_bucket_start != current_bucket_end) {
	    	    
	    // first relax higher degree vertices
	    for (typename Bucket::size_type i = current_bucket_start + tid ; 
		 i < current_bucket_end ; 
		 i+= nthreads) {
	      vertex_component_data vcd = (*buckets[current_bucket])[i];	    
	      Vertex v = vcd.first;
	      Vertex lv = logical_id(v);
	      Component dv = get(components, v);

	      if (dv < vcd.second)
		continue;
	      // we only need to relax if the component id
	      // has not changed from its logical id
	      bool haslowernbr = false;
	      std::set<Vertex> adjacencies;
	      BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
		Vertex u = target(e, g);
		Vertex lu = logical_id(target(e, g));
		if (lu != lv) { // ignore self-loops
		  if (lu > dv) {
		    if (adjacencies.insert(u).second) {
		      //relax_msg.send(vertex_component_data(u, dv));
#ifdef PRINT_DEBUG
		      std::cout << "  Relax " << get(get(vertex_local, g), v) << "@" 
				<< get(get(vertex_owner, g), v) << "->"
				<< get(get(vertex_local, g), u) << "@" 
				<< get(get(vertex_owner, g), u) << "  component = "
				<< dv << std::endl;
#endif
		    }
		  } else if (lu < dv) {
		    haslowernbr = true;
		    break;
		  }
		}
	      }

	      if (!haslowernbr) {
		typename std::set<Vertex>::iterator ite = adjacencies.begin();
		for(; ite != adjacencies.end(); ++ite) {
		  relax_msg.send(vertex_component_data((*ite), dv));
#ifdef PRINT_DEBUG
		  std::cout << "  Relax " << get(get(vertex_local, g), v) << "@" 
			    << get(get(vertex_owner, g), v) << "->"
			    << get(get(vertex_local, g), u) << "@" 
			    << get(get(vertex_owner, g), u) << "  component = "
			    << dv << std::endl;
#endif
		}
	      }


	    }
	  
	    // Wait for all threads to finish current bucket
	    t_bar->wait();
	    current_bucket_start = current_bucket_end;
	    current_bucket_end = buckets[current_bucket]->size();       
	    t_bar->wait();
	  
#ifdef PRINT_DEBUG
	    if (tid == 0)
	      std::cout << transport.rank() << ": end-processing bucket "
			<< current_bucket << " [" << current_bucket_start << ", " << current_bucket_end << ")\n";
#endif
	  } // end of (current_bucket_start != current_bucket_end) 
	} // No more re-insertions into this bucket to process (end of epoch)
      
	count_epoch += 1;
	// If all processes are done with the current bucket:
	//   1. clear the current bucket
	//   2. process low degree vertices
	//   3. find the next bucket to work on 
	if (all_starting_sizes == 0) {
	  
	  if (tid == 0) buckets[current_bucket]->clear();
	  t_bar->wait();


#ifdef PRINT_DEBUG
	  std::cout << tid << "@" << transport.rank() << ": Current bucket " << current_bucket
		    << " size " << buckets[current_bucket]->size() << std::endl;
#endif
	  assert(buckets[current_bucket]->size() == 0);


	  // find next bucket with work and update current_level 
	  BucketIndex old_bucket = current_bucket;
          
	  t_bar->wait();
	  if (tid == 0) {
	    find_next_bucket();

	    if (current_bucket != (std::numeric_limits<BucketIndex>::max)()) {
	      current_level += current_bucket - old_bucket;
	      ++buckets_processed;
	    }
	  }
	  t_bar->wait();
  
	  break; // exit the do..while loop for the current bucket
	} else // Update the end of the bucket and continue processing it
	  current_bucket_end = buckets[current_bucket]->size();
        
      } while(true);
      
      // If there are no non-empty buckets in any process, we're done
    } while(current_bucket != (std::numeric_limits<BucketIndex>::max)());
  }

  end_time = get_time();
}

template<DELTA_STEPPING_CC_PARMS>
void
DELTA_STEPPING_CC_TYPE::relax(const vertex_component_data& vcd,
			      int tid) {

  Vertex v = vcd.first;
  Component d = vcd.second;
  using boost::parallel::val_compare_and_swap;

  Component old_comp = get(components, v), last_old_comp;
  while (d < old_comp) {
    last_old_comp = old_comp;
    old_comp = val_compare_and_swap(&components[v], old_comp, d);
    if (last_old_comp == old_comp) {

#ifdef PBGL2_PRINT_WORK_STATS
      if (last_old_comp < logical_id(v))
        work_stats.increment_invalidated(tid);
      else
        work_stats.increment_useful(tid);
#endif

      // Insert vertex into new bucket, note we don't remove it from any other
      // buckets it might be in 
      BucketIndex new_index = 
        static_cast<BucketIndex>((d - (current_level * delta)) / delta);

      new_index = (current_bucket + new_index) % num_buckets;
      buckets[new_index]->push_back(vcd);

      // new_components now == get(components, v) so we exit automatically
      // but this saves us the conditional test
      return;
    }
  }
 
#ifdef PBGL2_PRINT_WORK_STATS
  work_stats.increment_rejected(tid);
#endif

}

template<DELTA_STEPPING_CC_PARMS>
void
DELTA_STEPPING_CC_TYPE::find_next_bucket()
{
  using boost::parallel::all_reduce;
  using boost::parallel::minimum;

  BucketIndex old_bucket = current_bucket;
  BucketIndex max_bucket = (std::numeric_limits<BucketIndex>::max)();

  current_bucket = (current_bucket + 1) % buckets.size();
  while (current_bucket != old_bucket && buckets[current_bucket]->empty())
    current_bucket = (current_bucket + 1) % buckets.size();
  
  if (current_bucket == old_bucket) 
    current_bucket = max_bucket;
  
  // If we wrapped, project index past end of buckets to use min()
  if (current_bucket < old_bucket) current_bucket += buckets.size();
  
  all_reduce<BucketIndex, minimum<BucketIndex> > r(transport, minimum<BucketIndex>());
  current_bucket = r(current_bucket);
  
  // Map index back into range of buckets
  if (current_bucket != max_bucket) 
    current_bucket %= buckets.size();
}

template<DELTA_STEPPING_CC_PARMS>
struct DELTA_STEPPING_CC_TYPE::vertex_component_handler {
  
  vertex_component_handler() : self(NULL) {}
  vertex_component_handler(delta_stepping_cc& self) : self(&self) {}
  
  void operator() (const vertex_component_data& data) const {
    int tid = amplusplus::detail::get_thread_id();
#ifdef PBGL2_PRINT_WORK_STATS
  self->work_stats.increment_edges(tid);
#endif
    if (data.second < get(self->components, data.first))
      self->relax(data, tid);
  }

protected:
  delta_stepping_cc* self;
};


} } } // end namespace boost::graph::distributed

#endif // BOOST_GRAPH_DELTA_STEPPING_CC_HPP
