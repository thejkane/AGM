#ifndef AMPP_PERF_COUNTER_DRIVER_INC
#define AMPP_PERF_COUNTER_DRIVER_INC

#include <am++/am++.hpp>
#include <am++/mpi_transport.hpp>
#include <am++/counter_coalesced_message_type.hpp>
#include <mpi.h>

//
// Performance counters
//
#ifdef AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS

std::vector<time_type> epoch_times;
// These should really be atomic
typedef unsigned long long flush_type;
typedef amplusplus::detail::atomic<flush_type> atomic_flush_type;
#define FLUSH_MPI_TYPE MPI_UNSIGNED_LONG
std::unique_ptr<atomic_flush_type[]> flushes;
flush_type flushes_size;
std::vector<flush_type> all_flushes, cumulative_flushes;
atomic_flush_type full{}, messages_received{};
flush_type all_full{}, cumulative_full{}, all_messages{}, cumulative_messages{};

void print_and_clear_epoch_times()
{
  if(_RANK == 0) {
    std::cout << "There were " << epoch_times.size() << " epochs." << std::endl;
    std::cout << "Epoch times ";
    BOOST_FOREACH(time_type t, epoch_times) {
      std::cout << print_time(t) << " ";
    }
    std::cout << "\n";
  }
  epoch_times.clear();
}

void clear_buffer_stats() {
  for(unsigned int i = 0; i < flushes_size; ++i)
    flushes[i] = 0;
  full = 0;
  messages_received = 0;
}

void clear_cumulative_buffer_stats() {
  for(unsigned int i = 0; i < flushes_size; ++i)
    cumulative_flushes[i] = 0;
  cumulative_full = 0;
  cumulative_messages = 0;
}

void print_buffer_stats() {
  unsigned long long sum = 0;
  unsigned long long number_of_buffers = 0;

  flush_type tmp_flushes[flushes_size];
  for(int i = 0; i < flushes_size; ++i)
    tmp_flushes[i] = flushes[i].load();
  flush_type tmp_full = full.load();
  flush_type tmp_messages = messages_received.load();

  MPI_Allreduce(&tmp_flushes, &all_flushes.front(), flushes_size, FLUSH_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&tmp_full, &all_full, 1, FLUSH_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&tmp_messages, &all_messages, 1, FLUSH_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD);

  if(_RANK == 0) std::cout << "Flushes: ";
  for(unsigned int i = 0; i < all_flushes.size(); ++i) {
    if(all_flushes[i] != 0) {
      if(_RANK == 0) std::cout << i+1 << "->" << all_flushes[i] << "  ";
      sum += (i + 1) * all_flushes[i];
      number_of_buffers += all_flushes[i];
      cumulative_flushes[i] += all_flushes[i];
    }
  }
  cumulative_full += all_full;
  cumulative_messages += all_messages;
  if(_RANK == 0) {
    std::cout << std::endl;
    std::cout << "There were " << number_of_buffers << " incomplete flushes with the average size of a flush of " << (number_of_buffers != 0 ? sum/number_of_buffers : 0) << std::endl;
    std::cout << "Full buffers: " << all_full << std::endl;
    std::cout << "Messages: " << all_messages << std::endl;
  }
}

namespace amplusplus {
  namespace performance_counters {
    void hook_flushed_message_size(amplusplus::rank_type dest, size_t count, size_t elt_size)
    {
      if(count <= flushes_size)
	flushes[count-1]++;

#ifdef PRINT_STATS
      time_type t = get_time();
      fprintf(stderr, "Flush: %d bytes to %d at %f\n", (count * elt_size), dest, t);
#endif
    }

    void hook_full_buffer_send(amplusplus::rank_type dest, size_t count, size_t elt_size)
    {
      full++;
#ifdef PRINT_STATS
      time_type t = get_time();
      fprintf(stderr, "Full buffer: %d bytes to %d at %f\n", (count * elt_size), dest, t);
#endif
    }

    void hook_message_received(amplusplus::rank_type src, size_t count, size_t elt_size)
    {
      messages_received += count;
#ifdef PRINT_STATS
      time_type t = get_time();
      fprintf(stderr, "%d: Message received: %d bytes from %d at %f\n", _RANK, (count * elt_size), src, t);
#endif
    }

    void hook_begin_epoch(amplusplus::transport& trans)
    {
      epoch_times.push_back(get_time());
    }

    void hook_epoch_finished(amplusplus::transport& trans)
    {
      epoch_times.back() = get_time() - epoch_times.back();
    }

  }
}
#endif // AMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS

// Hit rates
#ifdef AMPLUSPLUS_PRINT_HIT_RATES
unsigned long long cumulative_hits, cumulative_tests;

void print_and_accumulate_cache_stats(std::pair<unsigned long long, unsigned long long> stats) {
  unsigned long long hits, tests;
  MPI_Allreduce(&stats.first, &hits, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&stats.second, &tests, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  cumulative_hits += hits;
  cumulative_tests += tests;
  if(_RANK == 0) {
    std::cout << "Cache hit rates: " << hits << " hits, " << tests << " tests. Success ratio: " << (double)hits / (double)tests << "." << std::endl;
  }
}
#endif // AMPLUSPLUS_PRINT_HIT_RATES


#endif
