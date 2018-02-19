#ifndef __GIZMO_STAT__
#define __GIZMO_STAT__

#include <iostream>

namespace boost { namespace graph { namespace agm {
      
class stat_reader {

private:
  uint64_t useful_work;
  uint64_t rejected_work;
  uint64_t invalidated_work;

public:
  stat_reader() : useful_work(0), rejected_work(0), invalidated_work(0) {}

  void print() {
    std::cout << "[INFO] Useful work : " << useful_work << std::endl;
    std::cout << "[INFO] Rejeceted work : " << rejected_work << std::endl;
    std::cout << "[INFO] Invalidated work : " << invalidated_work << std::endl;
  }

  void increment_invalid() {
    ++invalidated_work;
  }

  void increment_reject() {
    ++rejected_work;
  }

  void increment_useful() {
    ++useful_work;
  }

  void reset() {
    useful_work = 0;
    rejected_work = 0;
    invalidated_work = 0;
  }
};

}}}
#endif
