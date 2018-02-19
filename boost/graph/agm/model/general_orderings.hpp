#ifndef __AGM_GENERAL_ORDERINGS__
#define __AGM_GENERAL_ORDERINGS__

#include <assert.h>
#include <tuple>
#include <iostream>
#include <map>

#define CHAOTIC_ORDERING_T boost::graph::agm::chaotic
#define DIJKSTRA_ORDERING_T boost::graph::agm::dijkstra<1>
#define DIJKSTRA_ORDERING_STD_PQ_T boost::graph::agm::dijkstra_std_pq<1>
#define DELTA_ORDERING_T boost::graph::agm::delta_ord<1, int>
#define LEVEL_ORDERING_T boost::graph::agm::level<1>
#define KLEVEL_ORDERING_T boost::graph::agm::klevel<1>
#define SSSP_KLEVEL_ORDERING_T boost::graph::agm::klevel<2>

namespace boost { namespace graph { namespace agm {

enum eagm_ordering {
  enum_chaotic,
  enum_dijkstra,
  enum_dijkstra_std_pq,
  enum_delta,
  enum_vertex_delta,
  enum_level,
  enum_klevel  
};
      
class base_ordering {
public:  
  virtual eagm_ordering name()=0;

  static std::string get_name(eagm_ordering e) {
    switch(e) {
    case eagm_ordering::enum_chaotic:
      return "chaotic";
    case eagm_ordering::enum_dijkstra:
      return "dijkstra";
    case eagm_ordering::enum_dijkstra_std_pq:
      return "dijkstra-std-pq";
    case eagm_ordering::enum_delta:
      return "delta";
    case eagm_ordering::enum_level:
      return "level";
    case eagm_ordering::enum_klevel:
      return "klevel";
    case eagm_ordering::enum_vertex_delta:
      return "vertex-delta";
    default:
      std::cout << "Invalid EAGM ordering. " << std::endl;
      break;
    }
  }

  static eagm_ordering get_ordering_by_name(std::string s) {
    if (s == "chaotic")
      return eagm_ordering::enum_chaotic;
    else if(s == "dijkstra")
      return eagm_ordering::enum_dijkstra;
    else if(s == "dijkstra-std-pq")
      return eagm_ordering::enum_dijkstra_std_pq;
    else if (s == "delta")
      return eagm_ordering::enum_delta;
    else if (s == "level")
      return eagm_ordering::enum_level;
    else if (s == "klevel")
      return eagm_ordering::enum_klevel;
    else if (s == "vertex-delta")
      return eagm_ordering::enum_vertex_delta;
    else {
      std::cout << "Invalid ordering name. Available : chaotic, dijkstra, dijkstra-std-pq, delta, level" << std::endl;
      assert(false);
    }
  }
};
      
// strict weak orderings
//================== Chaotic  =====================
struct chaotic : public base_ordering {
public:
  static const eagm_ordering ORDERING_NAME = eagm_ordering::enum_chaotic;
  
  template <typename T>
  bool operator()(T i, T j) {
    return false;
  }

  eagm_ordering name() {
    return ORDERING_NAME;
  }
};


//================== Chaotic  =====================

//================== Dijkstra  =====================
template<int index>
struct dijkstra : public base_ordering {
public:
  static const eagm_ordering ORDERING_NAME = eagm_ordering::enum_dijkstra;
  
  template <typename T>
  bool operator()(T i, T j) {
    return (std::get<index>(i) < std::get<index>(j));
  }

  eagm_ordering name() {
    return ORDERING_NAME;
  }
};

template<int index>
struct dijkstra_std_pq :public base_ordering {
public:
  static const eagm_ordering ORDERING_NAME = eagm_ordering::enum_dijkstra_std_pq;
  
  template <typename T>
  bool operator()(T i, T j) {
    return (std::get<index>(i) > std::get<index>(j));
  }

  eagm_ordering name() {
    return ORDERING_NAME;
  }
};
      

//================== Dijkstra  =====================

//================== Delta  =====================      
template<int index, typename delta_t>
struct delta_ord : public base_ordering {

private:
  delta_t delta;
  
public:
  static const eagm_ordering ORDERING_NAME = eagm_ordering::enum_delta;  
  delta_ord() : delta(10) {}
  delta_ord(delta_t _d) : delta(_d) {}
  delta_ord(const delta_ord& _do) : delta(_do.delta) {}
  
  template <typename T>
  bool operator()(T i, T j) {
    return ((std::get<index>(i)/delta) < (std::get<index>(j)/delta));
  }

  eagm_ordering name() {
    return ORDERING_NAME;
  }
};

template<typename delta_t,
         typename id_distribution>
struct cc_delta_ord : public base_ordering {

private:
  delta_t delta;
  id_distribution& idd;
  
public:
  static const eagm_ordering ORDERING_NAME = eagm_ordering::enum_vertex_delta;  
  cc_delta_ord(delta_t _d, id_distribution& _idd) : delta(_d),
                                                 idd(_idd){}
  cc_delta_ord(const cc_delta_ord& _do) : delta(_do.delta),
                                    idd(_do.idd){}
  
  template <typename T>
  bool operator()(T i, T j) {
    return ((idd(std::get<0>(i))/delta) < (idd(std::get<0>(j))/delta));
  }

  eagm_ordering name() {
    return ORDERING_NAME;
  }
};      
//================== Delta  =====================

//================== Level (same as Dijkstra)  =====================
template<int index>
struct level : public base_ordering {
public:
  static const eagm_ordering ORDERING_NAME = eagm_ordering::enum_level;    
  template <typename T>
  bool operator()(T i, T j) {
    return (std::get<index>(i) < std::get<index>(j));
  }

  eagm_ordering name() {
    return ORDERING_NAME;
  }
};


template<int index>
struct klevel : public base_ordering {
private:
  int k;
  
public:
  static const eagm_ordering ORDERING_NAME = eagm_ordering::enum_klevel;

  klevel():k(1){}
  klevel(int _k):k(_k){}  
  
  klevel(const klevel &o2) {k = o2.k;}
  
  template <typename T>
  bool operator()(T i, T j) {
    return ((std::get<index>(i)/k) < (std::get<index>(j)/k));
  }

  void set_value(int _k) {
    k = _k;
  }

  eagm_ordering name() {
    return ORDERING_NAME;
  }
};
      
class ordering_type_map {

public:

  typedef std::map<eagm_ordering, base_ordering*> ordering_map_t;
  
  static void init() {
    base_ordering* pord = new CHAOTIC_ORDERING_T;
    ordering_map.insert(std::make_pair(pord->name(), pord));
    pord = new DIJKSTRA_ORDERING_T;
    ordering_map.insert(std::make_pair(pord->name(), pord));
    pord = new DIJKSTRA_ORDERING_STD_PQ_T;
    ordering_map.insert(std::make_pair(pord->name(), pord));
    pord = new DELTA_ORDERING_T;
    ordering_map.insert(std::make_pair(pord->name(), pord));
    pord = new LEVEL_ORDERING_T;    
    ordering_map.insert(std::make_pair(pord->name(), pord));
    pord = new KLEVEL_ORDERING_T;    
    ordering_map.insert(std::make_pair(pord->name(), pord));    
    
  }

  static base_ordering* get_ordering(std::string _s) {
    auto find = ordering_map.find(base_ordering::get_ordering_by_name(_s));
    if (find != ordering_map.end()) {
      return find->second;
    } else
      return NULL;
  }

  static ordering_map_t create_ordering_map() {
    ordering_map_t map;
    return map;
  }

  static void print() {
    for (auto it = ordering_map.begin(); it != ordering_map.end(); ++it) {
      std::cout << it->first << std::endl;
    }
  }

  static void clear() {
    for (auto it = ordering_map.begin(); it != ordering_map.end(); ++it) {
      delete it->second;
    }

    ordering_map.clear();
  }

private:
  static ordering_map_t ordering_map;
  
};

}}}
#endif
