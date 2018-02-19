#ifndef __EAGM_CONFIG_H__
#define __EAGM_CONFIG_H__

namespace boost { namespace graph { namespace agm {
//=====================================================//
// EAGM Configuration for Static Spatial Structure (default)
//=====================================================//


struct buffer_container {};
struct pq_container {};
      
template<typename T1,
         typename T2,
         typename T3,
         typename T4,
         typename node_container=buffer_container,
         typename numa_container=buffer_container>
class eagm_configs {
  
public:
  typedef T1 global_ordering_t;
  typedef T2 node_ordering_t;
  typedef T3 numa_ordering_t;
  typedef T4 thread_ordering_t;
  typedef node_container node_container_t;
  typedef numa_container numa_container_t;  

  eagm_configs(T1 _go,
               T2 _no,
               T3 _nuo,
               T4 _to,
	       bool _forward_scheduler=true) : global_ord(_go),
					       node_ord(_no),
					       numa_ord(_nuo),
					       thread_ord(_to),
					       forward_scheduler(_forward_scheduler){}

  eagm_configs(const eagm_configs& _ec): global_ord(_ec.global_ord),
                                         node_ord(_ec.node_ord),
                                         numa_ord(_ec.numa_ord),
                                         thread_ord(_ec.thread_ord),
					 forward_scheduler(_ec.forward_scheduler){}
  

  void print() {
    std::cout << "[Global:" << base_ordering::get_name(global_ord.name()) << "]"
              << "-->"
              << "[Node:" << base_ordering::get_name(node_ord.name()) << "]"
              << "-->"
              <<"[Numa:" << base_ordering::get_name(numa_ord.name()) << "]"
              << "-->"
              <<"[Thread:" << base_ordering::get_name(thread_ord.name()) << "]"
              << std::endl;
  }
  
public:
  global_ordering_t global_ord;
  node_ordering_t node_ord;
  numa_ordering_t numa_ord;
  thread_ordering_t thread_ord;
  bool forward_scheduler;
};

template<typename GT1,
         typename GT2,
         typename GT3,
         typename GT4>
static auto create_eagm_config(GT1 _t1, GT2 _t2, GT3 _t3, GT4 _t4) -> decltype(eagm_configs<GT1, GT2, GT3, GT4>(_t1, _t2, _t3, _t4)) {
    return eagm_configs<GT1, GT2, GT3, GT4>(_t1, _t2, _t3, _t4);    
  }      
      
enum spatial_level {
  global = 0,
  node = 1,
  numa = 2,
  thread = 3,
  nil = 4
};

      
//=====================================================//
// EAGM Configuration for Dynamic Spatial Structure
//=====================================================//      

      
template<typename T1,
         typename T2,
         typename T3,
         typename T4>
class dynmic_eagm_configs {
public:
  typedef T1 global_ordering_t;
  typedef T2 node_ordering_t;
  typedef T3 numa_ordering_t;
  typedef T4 thread_ordering_t;

  dynmic_eagm_configs(T1 _go,
               T2 _no,
               T3 _nuo,
               T4 _to) : global_ord(_go),
                         node_ord(_no),
                         numa_ord(_nuo),
                         thread_ord(_to){}

  dynmic_eagm_configs(const dynmic_eagm_configs& _ec): global_ord(_ec.global_ord),
                                         node_ord(_ec.node_ord),
                                         numa_ord(_ec.numa_ord),
                                         thread_ord(_ec.thread_ord),
                                         ranks(_ec.ranks),
                                         threads(_ec.threads),
                                         numa_domains(_ec.numa_domains){}

public:
  global_ordering_t global_ord;
  node_ordering_t node_ord;
  numa_ordering_t numa_ord;
  thread_ordering_t thread_ord;
  int ranks;
  int threads;
  int numa_domains;

  spatial_level get_child_spatial_level(spatial_level parent) {
    return optimized_spatial_hierarchy[parent];
  }


  void print() {
    std::cout << "ranks = " << ranks << std::endl;
    std::cout << "threads = " << threads << std::endl;
    std::cout << "numa_domains = " << numa_domains << std::endl;

    std::cout << "printing the hierarchy ..." << std::endl;

    std::cout << "Global --> ";

    int nextlevel = global;

    while(optimized_spatial_hierarchy[nextlevel] != nil) {
      if (optimized_spatial_hierarchy[nextlevel] == thread) {
        std::cout << " Thread --> ";
        nextlevel = thread;
      } else if (optimized_spatial_hierarchy[nextlevel] == node) {
        std::cout << " Node --> ";
        nextlevel = node;
      } else if (optimized_spatial_hierarchy[nextlevel] == numa) {
        std::cout << " Numa --> ";
        nextlevel = numa;
      }
    }

    std::cout << " nil" << std::endl;
  }
  
  void optimize() {
    
    optimized_spatial_hierarchy.resize(4, nil);
    optimized_spatial_hierarchy.reserve(4);
    
    // condition 1 : are all orderings chaotic ? map global to nil
    bool chaotic_node = std::is_same<node_ordering_t,
                                     CHAOTIC_ORDERING_T>::value;
    bool chaotic_numa = std::is_same<numa_ordering_t,
                                     CHAOTIC_ORDERING_T>::value;
    bool chaotic_thread = std::is_same<thread_ordering_t,
                                       CHAOTIC_ORDERING_T>::value;

    if (chaotic_node) {
      if (chaotic_numa) {
        if (chaotic_thread) {
          // all are chaotic
          optimized_spatial_hierarchy[global] = nil;
        } else {
          // thread is not chaotic but node and numa is
          optimized_spatial_hierarchy[global] = thread;
        }
      } else {
        // node is chaotic, numa is not
        optimized_spatial_hierarchy[global] = numa;

        if (chaotic_thread)
          optimized_spatial_hierarchy[numa] = nil;
        else
          optimized_spatial_hierarchy[numa] = thread;
      }
    } else {
      // node is not chaotic
      optimized_spatial_hierarchy[global] = node;

      if (chaotic_numa) {
        // numa is chaotic
        if (chaotic_thread) {
          // numa chaotic thread chaotic
          optimized_spatial_hierarchy[node] = nil;
        } else {
          // numa chaotic and thread not chaotic
          optimized_spatial_hierarchy[node] = thread;
        }
      } else {
        // numa is not chaotic
        optimized_spatial_hierarchy[node] = numa;

        if (chaotic_thread)
          optimized_spatial_hierarchy[numa] = nil;
        else
          optimized_spatial_hierarchy[numa] = thread;
      }
    }    
  }

private:
  std::vector<spatial_level> optimized_spatial_hierarchy;
};

}}}
#endif
