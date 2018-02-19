#ifndef __AGM_BFS_HPP
#define __AGM_BFS_HPP

#include "../models/algorithm/agm.hpp"
#include "../models/machine/access_metadata.hpp"
#include "../models/machine/shared_memory.hpp"
#include "../stats/stat.hpp"

#include <boost/config.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/iteration_macros.hpp>
//#include <boost/graph/dijkstra_shortest_paths.hpp>

#include <boost/graph/iteration_macros.hpp>

// definition of the BFS work item set
typedef int Level;

template<typename Graph>
class bfs_simulator {

  typedef typename boost::graph_traits < Graph >::vertex_descriptor Vertex;

  // work item set definition
  // Every work item must have the time
  // destination, level, timestamp, source
  typedef std::tuple<Vertex, Level, uint64_t, Vertex> WorkItem;

  // the processing function  
  template<typename State>
  struct bfs_pf {

  public:
    bfs_pf(Graph& _rg, State& _st, stat_reader& _sr) : g(_rg),
						       vlevel(_st),
						       stats(_sr){}

  private:
    Graph& g;
    State& vlevel;
    stat_reader& stats;

  public:
    void operator()(WorkItem& wi,
		    std::vector<WorkItem>& outset,
		    uint64_t& ts) {

      auto new_distance = std::get<1>(wi);
      std::string color = "\033["+std::to_string(0)+"m";
      if (new_distance > 0) {
	color = "\033[0;"+std::to_string(new_distance+29)+"m";
      }

      //#ifdef PRINT_DEBUG
      std::cout << color << "[(" << std::get<0>(wi) << "," << std::get<3>(wi) << ")" << std::get<1>(wi) << "] ";
      //      std::cout << "processing v : [(" << std::get<0>(wi) << ", " << std::get<3>(wi) << ")" << std::get<1>(wi) << ", " << std::get<2>(wi) << ")" << std::endl;
      //#endif
      // accessed a state
      //util::increment_time_no_cb(wi);
      ++ts;
      if (std::get<1>(wi) < vlevel[std::get<0>(wi)]) {

	if ((vlevel[std::get<0>(wi)]) == INT_MAX) {
	  debug("=========================== USEFUL ==============================");
	  stats.increment_useful();
	} else
	  stats.increment_invalid();

	// updating a state
	//	util::increment_time_no_cb(wi);
	++ts;
        vlevel[std::get<0>(wi)] = std::get<1>(wi);

	// reading the index array
	//	util::increment_time_no_cb(wi);
	++ts;
        BGL_FORALL_OUTEDGES_T(std::get<0>(wi), e, g, Graph) {
	  // for each adjacent increment
	  //	  util::increment_time_no_cb(wi);
	  ++ts;
          Vertex v = boost::target(e, g);
	  
	  // std::cout << "Time before creation : " << std::get<2>(wi) << std::endl;
	  WorkItem generated(v, (std::get<1>(wi)+1), std::get<2>(wi), std::get<0>(wi));
	  // std::cout << "Time after creation : " << std::get<2>(generated) << std::endl;
	  outset.push_back(generated);
        }
      } else {
	debug("Rejected work ...");
	stats.increment_reject();
      }
    }
  };

  // states property map

public:
  template<typename DistMap>
  bool verify(Graph& g, Vertex s, DistMap& dists) {
    info("Verifying distances ....");

    BGL_FORALL_VERTICES_T(v, g, Graph) {
      BGL_FORALL_OUTEDGES_T(v, e, g, Graph) {
	if (dists[source(e, g)] == INT_MAX) {
	  if (dists[target(e, g)] != INT_MAX) {
	    std::cout << "[ERROR] Source unvisited, therefore target must also be unvisited. But Target -- (" << target(e, g) << ", " << dists[target(e, g)] << ") -> Source -- ("
		      << source(e, g) << ", " << dists[source(e, g)] << ") " << std::endl;
	    std::abort();
	  }
	} else if (dists[target(e, g)] > (dists[source(e, g)] + 1)) {
	  std::cout << "[ERROR] Target distance is greater than source distance + 1: Target -- (" << target(e, g) << ", " << dists[target(e, g)] << ") -> Source -- ("
		    << source(e, g) << ", " << dists[source(e, g)] << ") " << std::endl;	  
	  std::abort();
	}
      }
    }

    info("Verification successful ...");
    return true;
  }

  template<typename Ordering,
	   typename MachineModelGen>
  void simulate(Graph& g, 
		stat_reader& sr, 
		gizmo_config& config,
		bool _verify=true) {

    info("Inside BFS simulator ...");
    info("Creating the state ...");
    // State
    typedef typename boost::property_map<Graph, boost::vertex_index_t>::type VertexIndexMap;
    std::vector<int> distmap(num_vertices(g), INT_MAX);
    typedef boost::iterator_property_map<typename std::vector<int>::iterator, VertexIndexMap> DistMap;
    DistMap distance_state(distmap.begin(), get(boost::vertex_index, g));

    // source
    Vertex source = 2;

    info("Setting the initial work item set ...");
    // Initial work item set
    std::vector<WorkItem> initial;
    initial.push_back(WorkItem(source, 0, 0, source));

    info("Setting the ordering ...");
    // Ordering
    typedef typename Ordering::strict_ordering StrictWeakOrdering;

    info("Setting the processing function ...");
    // Processing funcion
    typedef bfs_pf<DistMap> ProcessingFunction;
    ProcessingFunction pf(g, distance_state, sr);

    info("Setting the machine model ...");
    // Machine model 
    typedef typename MachineModelGen::template inner<WorkItem, StrictWeakOrdering, 
						     ProcessingFunction>::type MachineModel;
    info("Instantiating machine model ...");
    MachineModel model(boost::num_vertices(g), pf);
    info("Initializing the machine model ...");
    model.initialize(config);

    info("Setting the callback ...");
    // set callback -- this is a must
    util::set_callback(&model);

    // BFS algorithm
    typedef agm<WorkItem, StrictWeakOrdering, MachineModel> bfs_chaotic_agm_t;

    bfs_chaotic_agm_t bfschaotic;
    
    info("Invoking BFS simulator with chaotic ordering ...");
    bfschaotic(model, initial);

    if (_verify) {
      verify(g, source, distance_state);
    }

    std::cout << "[INFO] Execution time : " << model.execution_time() << std::endl;
    std::cout << "[INFO] Fastest worker time : " << model.fastest_execution_time() << std::endl;
  }
};


#endif
