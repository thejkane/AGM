#ifndef __EAGM_BUCKET_TRAITS__
#define __EAGM_BUCKET_TRAITS__

#include <boost/thread/locks.hpp>
#include <boost/parallel/append_buffer.hpp>
#include <boost/graph/agm/model/general_orderings.hpp>
#include <boost/graph/agm/util/eagm_config.hpp>

namespace boost { namespace graph { namespace agm {

struct spatial_global_tag {};
struct spatial_node_tag {};
struct spatial_numa_tag {};
struct spatial_thread_tag {};
struct spatial_global_buffer_tag {};      
struct spatial_node_buffer_tag {};
// We need this tag to handle concurrent priority queue
// and append buffers
struct spatial_node_select_pq_or_buffer_tag {};
// We need this tag to handle concurrent priority queue
// and append buffers for numa domains.
struct spatial_numa_select_pq_or_buffer_tag {};            
struct spatial_numa_buffer_tag {};
struct spatial_global_all_chaotic_buffer_tag {};                  

struct spatial_levels_one {}; // e.g., Global--> Nil
struct spatial_levels_two {}; // e.g., Global --> Thread --> Nil
struct spatial_levels_three {}; // e.g., Global --> Node --> Numa
struct spatial_levels_four {}; // e.g., Global --> Node --> Numa --> Thread

struct BufferOrdering {};
struct leaf_level_eagm_trait {};      

// EAGM traits
template<typename Order,
         typename SpatialTag,
         typename EAGMConfig,
         typename Runtime,
         typename NextLevelEagmTrait=boost::graph::agm::leaf_level_eagm_trait>
struct eagm_bucket_traits {
  typedef Order order_t;
  typedef SpatialTag spatial_t;
  typedef EAGMConfig eagm_config;
  typedef Runtime runtime;
  typedef NextLevelEagmTrait next_level_trait_t;
};
            
// Global --> Node --> Numa --> Thread --> nil
template<typename GlobalOrdering,
         typename NodeOrdering,
         typename NumaOrdering,
         typename ThreadOrdering,
         typename EAGMConfig,
         typename Runtime>
struct eagm_traits {

  static_assert(std::is_same<GlobalOrdering, typename EAGMConfig::global_ordering_t>::value,
                "EAGM configuration global ordering and trait global ordering must be the same.");
  static_assert(std::is_same<NodeOrdering, typename EAGMConfig::node_ordering_t>::value,
                "EAGM configuration node ordering and trait node ordering must be the same.");
  static_assert(std::is_same<NumaOrdering, typename EAGMConfig::numa_ordering_t>::value,
                "EAGM configuration numa ordering and trait numa ordering must be the same.");
  static_assert(std::is_same<ThreadOrdering, typename EAGMConfig::thread_ordering_t>::value,
                "EAGM configuration thread ordering and trait thread ordering must be the same.");
  
  // Thread --> Nil
  typedef eagm_bucket_traits<ThreadOrdering,
                             spatial_thread_tag,
                             EAGMConfig,
                             Runtime> third_level_bucket_trait;

  // Numa --> Thread
  typedef eagm_bucket_traits<NumaOrdering,
                             spatial_numa_tag,
                             EAGMConfig,
                             Runtime,
                             third_level_bucket_trait> second_level_bucket_trait;

  // Node --> Numa
  typedef eagm_bucket_traits<NodeOrdering,
                             spatial_node_tag,
                             EAGMConfig,
                             Runtime,
                             second_level_bucket_trait> first_level_bucket_trait;
  
  // Global --> Node
  typedef eagm_bucket_traits<GlobalOrdering,
                             spatial_global_tag,
                             EAGMConfig,
                             Runtime,
                             first_level_bucket_trait> root_level_bucket_trait;

};

// Global --> Node --> Numa --> Nil
template<typename GlobalOrdering,
         typename NodeOrdering,
         typename NumaOrdering,
         typename EAGMConfig,
         typename Runtime>
struct eagm_traits<GlobalOrdering,
                   NodeOrdering,
                   NumaOrdering,
                   CHAOTIC_ORDERING_T,
                   EAGMConfig,
                   Runtime> {

  static_assert(std::is_same<GlobalOrdering, typename EAGMConfig::global_ordering_t>::value,
                "EAGM configuration global ordering and trait global ordering must be the same.");
  static_assert(std::is_same<NodeOrdering, typename EAGMConfig::node_ordering_t>::value,
                "EAGM configuration node ordering and trait node ordering must be the same.");
  static_assert(std::is_same<NumaOrdering, typename EAGMConfig::numa_ordering_t>::value,
                "EAGM configuration numa ordering and trait numa ordering must be the same.");


  // Note that we are adding a new level
  // buffer --> nil
  typedef eagm_bucket_traits<BufferOrdering,
                             spatial_numa_buffer_tag,
                             EAGMConfig,
                             Runtime> buffer_level_bucket_trait;

  
  // Numa --> Buffer
  typedef eagm_bucket_traits<NumaOrdering,
                             spatial_numa_select_pq_or_buffer_tag,
                             EAGMConfig,
                             Runtime,
                             buffer_level_bucket_trait> second_level_bucket_trait;

  // Node --> Numa
  typedef eagm_bucket_traits<NodeOrdering,
                             spatial_node_tag,
                             EAGMConfig,
                             Runtime,
                             second_level_bucket_trait> first_level_bucket_trait;
  
  // Global --> Node
  typedef eagm_bucket_traits<GlobalOrdering,
                             spatial_global_tag,
                             EAGMConfig,
                             Runtime,
                             first_level_bucket_trait> root_level_bucket_trait;

};

// Global-->Node-->Thread-->Nil
template<typename GlobalOrdering,
         typename NodeOrdering,
         typename ThreadOrdering,
         typename EAGMConfig,
         typename Runtime>
struct eagm_traits<GlobalOrdering,
                   NodeOrdering,
                   CHAOTIC_ORDERING_T,
                   ThreadOrdering,
                   EAGMConfig,
                   Runtime>{

  static_assert(std::is_same<GlobalOrdering, typename EAGMConfig::global_ordering_t>::value,
                "EAGM configuration global ordering and trait global ordering must be the same.");
  static_assert(std::is_same<NodeOrdering, typename EAGMConfig::node_ordering_t>::value,
                "EAGM configuration node ordering and trait node ordering must be the same.");
  static_assert(std::is_same<ThreadOrdering, typename EAGMConfig::thread_ordering_t>::value,
                "EAGM configuration thread ordering and trait thread ordering must be the same.");
  
  
  // Thread --> Nil
  typedef eagm_bucket_traits<ThreadOrdering,
                             spatial_thread_tag,
                             EAGMConfig,
                             Runtime> second_level_bucket_trait;

  // Node --> Thread
  typedef eagm_bucket_traits<NodeOrdering,
                             spatial_node_tag,
                             EAGMConfig,
                             Runtime,
                             second_level_bucket_trait> first_level_bucket_trait;
  
  // represents Global --> Node
  typedef eagm_bucket_traits<GlobalOrdering,
                             spatial_global_tag,
                             EAGMConfig,
                             Runtime,
                             first_level_bucket_trait> root_level_bucket_trait;  

};      

// Global-->Numa-->Thread-->Nil
template<typename GlobalOrdering,
         typename NumaOrdering,
         typename ThreadOrdering,
         typename EAGMConfig,
         typename Runtime>
struct eagm_traits<GlobalOrdering,
                   CHAOTIC_ORDERING_T,
                   NumaOrdering,
                   ThreadOrdering,
                   EAGMConfig,
                   Runtime>{

  static_assert(std::is_same<GlobalOrdering, typename EAGMConfig::global_ordering_t>::value,
                "EAGM configuration global ordering and trait global ordering must be the same.");
  static_assert(std::is_same<NumaOrdering, typename EAGMConfig::numa_ordering_t>::value,
                "EAGM configuration node ordering and trait node ordering must be the same.");
  static_assert(std::is_same<ThreadOrdering, typename EAGMConfig::thread_ordering_t>::value,
                "EAGM configuration thread ordering and trait thread ordering must be the same.");
  
  
  // Thread --> Nil
  typedef eagm_bucket_traits<ThreadOrdering,
                             spatial_thread_tag,
                             EAGMConfig,
                             Runtime> second_level_bucket_trait;

  // Numa --> Thread
  typedef eagm_bucket_traits<NumaOrdering,
                             spatial_numa_tag,
                             EAGMConfig,
                             Runtime,
                             second_level_bucket_trait> first_level_bucket_trait;
  
  // represents Global --> Node
  typedef eagm_bucket_traits<GlobalOrdering,
                             spatial_global_tag,
                             EAGMConfig,
                             Runtime,
                             first_level_bucket_trait> root_level_bucket_trait;  

};      

// Global-->Thread-->Nil
template<typename GlobalOrdering,
         typename ThreadOrdering,
         typename EAGMConfig,
         typename Runtime>
struct eagm_traits<GlobalOrdering,
                   CHAOTIC_ORDERING_T,
                   CHAOTIC_ORDERING_T,
                   ThreadOrdering,
                   EAGMConfig,
                   Runtime>{
  
  static_assert(std::is_same<GlobalOrdering, typename EAGMConfig::global_ordering_t>::value,
                "EAGM configuration global ordering and trait global ordering must be the same.");
  static_assert(std::is_same<ThreadOrdering, typename EAGMConfig::thread_ordering_t>::value,
                "EAGM configuration thread ordering and trait thread ordering must be the same.");
  
  // represents Thread-->Nil
  typedef eagm_bucket_traits<ThreadOrdering,
                             spatial_thread_tag,
                             EAGMConfig,
                             Runtime> first_level_bucket_trait;

  // represents Global --> Thread
  typedef eagm_bucket_traits<GlobalOrdering,
                             spatial_global_tag,
                             EAGMConfig,
                             Runtime,
                             first_level_bucket_trait> root_level_bucket_trait;  
};

// Global-->Node-->Nil
template<typename GlobalOrdering,
         typename NodeOrdering,
         typename EAGMConfig,
         typename Runtime>
struct eagm_traits<GlobalOrdering,
                   NodeOrdering,
                   CHAOTIC_ORDERING_T,
                   CHAOTIC_ORDERING_T,
                   EAGMConfig,
                   Runtime>{

  static_assert(std::is_same<GlobalOrdering, typename EAGMConfig::global_ordering_t>::value,
                "EAGM configuration global ordering and trait global ordering must be the same.");
  static_assert(std::is_same<NodeOrdering, typename EAGMConfig::node_ordering_t>::value,
                "EAGM configuration node ordering and trait node ordering must be the same.");

  // Note that we are adding a new level
  // buffer --> nil
  typedef eagm_bucket_traits<BufferOrdering,
                             spatial_node_buffer_tag,
                             EAGMConfig,
                             Runtime> buffer_level_bucket_trait;

  // node --> buffer
  typedef eagm_bucket_traits<NodeOrdering,
                             spatial_node_select_pq_or_buffer_tag, // Special tag
                             EAGMConfig,
                             Runtime,
                             buffer_level_bucket_trait> first_level_bucket_trait;
  
  // represents Global --> Node
  typedef eagm_bucket_traits<GlobalOrdering,
                             spatial_global_tag,
                             EAGMConfig,
                             Runtime,
                             first_level_bucket_trait> root_level_bucket_trait;
};

// Global-->Numa-->Nil
template<typename GlobalOrdering,
         typename NumaOrdering,
         typename EAGMConfig,
         typename Runtime>
struct eagm_traits<GlobalOrdering,
                   CHAOTIC_ORDERING_T,
                   NumaOrdering,
                   CHAOTIC_ORDERING_T,
                   EAGMConfig,
                   Runtime>{

  static_assert(std::is_same<GlobalOrdering, typename EAGMConfig::global_ordering_t>::value,
                "EAGM configuration global ordering and trait global ordering must be the same.");
  static_assert(std::is_same<NumaOrdering, typename EAGMConfig::numa_ordering_t>::value,
                "EAGM configuration numa ordering and trait numa ordering must be the same.");

  // buffer --> nil
  typedef eagm_bucket_traits<BufferOrdering,
                             spatial_numa_buffer_tag,
                             EAGMConfig,
                             Runtime> buffer_level_bucket_trait;

  // numa --> buffer
  typedef eagm_bucket_traits<NumaOrdering,
                             spatial_numa_select_pq_or_buffer_tag,
                             EAGMConfig,
                             Runtime,
                             buffer_level_bucket_trait> first_level_bucket_trait;
  // represents Global --> Numa
  typedef eagm_bucket_traits<GlobalOrdering,
                             spatial_global_tag,
                             EAGMConfig,
                             Runtime,
                             first_level_bucket_trait> root_level_bucket_trait;  
};

// Global-->Nil (the AGM)
template<typename GlobalOrdering,
         typename EAGMConfig,
         typename Runtime>
struct eagm_traits<GlobalOrdering,
                   CHAOTIC_ORDERING_T,
                   CHAOTIC_ORDERING_T,
                   CHAOTIC_ORDERING_T,
                   EAGMConfig,
                   Runtime>{

  static_assert(std::is_same<GlobalOrdering, typename EAGMConfig::global_ordering_t>::value,
                "EAGM configuration global ordering and trait global ordering must be the same.");

  // Note that we are adding a new level
  // buffer --> nil
  typedef eagm_bucket_traits<BufferOrdering,
                             spatial_global_buffer_tag,
                             EAGMConfig,
                             Runtime> buffer_level_bucket_trait;

  // Global --> Buffer
  typedef eagm_bucket_traits<GlobalOrdering,
                             spatial_global_tag,
                             EAGMConfig,
                             Runtime,
                             buffer_level_bucket_trait> root_level_bucket_trait;
};
      
// Pure Chaotic
// Global --> Buffer --> Nil
template<typename EAGMConfig,
         typename Runtime>
struct eagm_traits<CHAOTIC_ORDERING_T,
                   CHAOTIC_ORDERING_T,
                   CHAOTIC_ORDERING_T,
                   CHAOTIC_ORDERING_T,
                   EAGMConfig,
                   Runtime> {

  static_assert(std::is_same<CHAOTIC_ORDERING_T, typename EAGMConfig::global_ordering_t>::value,
                "EAGM configuration global ordering and trait global ordering must be the same.");
  static_assert(std::is_same<CHAOTIC_ORDERING_T, typename EAGMConfig::node_ordering_t>::value,
                "EAGM configuration node ordering and trait node ordering must be the same.");
  static_assert(std::is_same<CHAOTIC_ORDERING_T, typename EAGMConfig::numa_ordering_t>::value,
                "EAGM configuration numa ordering and trait numa ordering must be the same.");
  static_assert(std::is_same<CHAOTIC_ORDERING_T, typename EAGMConfig::thread_ordering_t>::value,
                "EAGM configuration thread ordering and trait thread ordering must be the same.");

    // buffer --> nil
  typedef eagm_bucket_traits<BufferOrdering,
                             spatial_global_all_chaotic_buffer_tag,
                             EAGMConfig,
                             Runtime> buffer_level_bucket_trait;

  
  // Global (Chaotic) --> Nil
  typedef eagm_bucket_traits<CHAOTIC_ORDERING_T,
                             spatial_global_tag,
                             EAGMConfig,
                             Runtime,
                             buffer_level_bucket_trait> root_level_bucket_trait;

};

}}}
#endif
