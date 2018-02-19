#ifndef __AGM_RUNTIME_COMMON__
#define __AGM_RUNTIME_COMMON__

#include <boost/property_map/property_map.hpp>
#include <boost/config.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/graph_traits.hpp>

template <typename OwnerMap, typename work_item>
struct work_item_owner {
  explicit work_item_owner(const OwnerMap& owner) : owner(owner) {}

  const OwnerMap& owner;
};

template <typename OwnerMap, typename work_item>
typename boost::property_traits<OwnerMap>::value_type
get(const work_item_owner<OwnerMap, work_item>& o, const work_item& data)
{ return get(o.owner, std::get<0>(data)); }

template<typename runtime>
class runtime_epoch: boost::noncopyable {
  runtime& rt;
  public:
  runtime_epoch(runtime& _rt): rt(_rt) {rt.begin_epoch(0);}
  ~runtime_epoch() {rt.end_epoch();}
};

template<typename runtime>
class runtime_epoch_value: boost::noncopyable {
  runtime& rt;
  const unsigned long& read_value;
  unsigned long& sum;
  public:
  runtime_epoch_value(runtime& _rt,
                     const unsigned long& _read_value,
                     unsigned long& sum)
    : rt(_rt), read_value(_read_value), sum(sum) {rt.begin_epoch(0);}
  ~runtime_epoch_value() {sum = rt.end_epoch_with_value(read_value);}
};

#endif
