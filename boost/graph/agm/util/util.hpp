#ifndef __AGM_UTILITIES__
#define __AGM_UTILITIES__
#include <iostream>
#include <tuple>
#include <limits>

template<typename T>
T get_init_val(T t) {
  return std::numeric_limits<T>::max();
}

bool get_init_val(bool t) {
  return false;
}

template<std::size_t I = 0, typename... Tp>
inline typename std::enable_if<I == sizeof...(Tp), void>::type
init_wi(std::tuple<Tp...>& t)
{}

template<std::size_t I = 0, typename... Tp>
inline typename std::enable_if<I < sizeof...(Tp), void>::type
init_wi(std::tuple<Tp...>& t) {
  using type = typename std::tuple_element<I, std::tuple<Tp...> >::type;
  type dummy;
  std::get<I>(t) = get_init_val(dummy);
  init_wi<I + 1, Tp...>(t);
}

template<std::size_t I = 0, typename... Tp>
inline typename std::enable_if<I == sizeof...(Tp), void>::type
iprint_wi(const std::tuple<Tp...>& t) {
}

template<std::size_t I = 0, typename... Tp>
inline typename std::enable_if<I < sizeof...(Tp), void>::type
iprint_wi(const std::tuple<Tp...>& t) {
    std::cout << std::get<I>(t) << ", ";
    iprint_wi<I + 1, Tp...>(t);
}

template<typename work_item_t>
void print_wi(const work_item_t& w) {
  std::cout << "(";
  iprint_wi(w);
  std::cout << ")";
}

#endif
