# Copyright (C) 2013 The Trustees of Indiana University
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met: 

# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer. 
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution. 

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# Authors: Marcin Zalewski

ACLOCAL_AMFLAGS=-I m4

nobase_include_HEADERS = \
	am++/am++.hpp \
	am++/basic_coalesced_message_type.hpp \
	am++/counter_coalesced_message_type.hpp \
	am++/size_coalesced_message_type.hpp \
	am++/lock_free_coalesced_message_type.hpp \
	am++/lock_free_coalesced_message_type_packed.hpp \
	am++/make_mpi_datatype.hpp \
	am++/mattern_channel_counting_termination_detector.hpp \
	am++/message_queue.hpp \
	am++/message_type_generators.hpp \
	am++/mpi_sinha_kale_ramkumar_termination_detector.hpp \
	am++/mpi_sinha_kale_ramkumar_termination_detector_bgp.hpp \
	am++/mpi_transport.hpp \
	am++/object_based_addressing.hpp \
	am++/one_sided.hpp \
	am++/performance_counters.hpp \
	am++/reductions.hpp \
	am++/scoped_epoch.hpp \
	am++/termination_detector.hpp \
	am++/traits.hpp \
	am++/transport.hpp \
	am++/detail/amortized_hook.hpp \
	am++/detail/append_buffer.hpp \
	am++/detail/buffer_cache.hpp \
	am++/detail/factory_wrapper.hpp \
	am++/detail/id_assigner.hpp \
	am++/detail/mpi_global_lock.hpp \
	am++/detail/mpi_pool.hpp \
	am++/detail/mpi_request_manager.hpp \
	am++/detail/nbc_check.hpp \
	am++/detail/signal.hpp \
	am++/detail/term_detect_level_manager.hpp \
	am++/detail/thread_support.hpp \
	am++/detail/type_info_map.hpp \
	am++/detail/typed_in_place_factory_owning.hpp \
	am++/detail/vector_of_noncopyable.hpp \
	am++/dummy_buffer_sorter.hpp

sources = \
	src/mpi_make_mpi_datatype.cpp \
	src/mpi_sinha_kale_ramkumar_termination_detector.cpp \
	src/mpi_sinha_kale_ramkumar_termination_detector_bgp.cpp \
	src/mpi_transport.cpp \
	src/termination_detector.cpp \
	src/thread_support.cpp \
	src/transport.cpp

lib_LTLIBRARIES = libampp.la
libampp_la_SOURCES = $(sources)
libampp_la_LDFLAGS = -version-info 0:0:0



check_HEADERS = \
	tests/distributed_erdos_renyi_generator.hpp \
	tests/rmat_graph_generator_faster.hpp \
	tests/splittable_ecuyer1988.hpp
check_PROGRAMS = \
	tests/test_bfs \
	tests/test_bfs_threaded \
	tests/test_matvec \
	tests/test_matvec_t \
	tests/test_pingpong \
	tests/test_ring  \
	tests/test_message_priority \
	tests/test_buffer

LDADD = libampp.la $(BOOST_LDFLAGS) $(BOOST_THREAD_LIB) $(BOOST_SYSTEM_LIB)
TESTS = \
        tests/test_bfs_threaded.exec1 \
        tests/test_bfs_threaded.exec2 \
	tests/test_bfs.exec \
	tests/test_matvec.exec \
	tests/test_matvec_t.exec \
	tests/test_pingpong.exec \
	tests/test_ring.exec \
	tests/test_message_priority.exec \
	tests/test_buffer.exec

AM_DEFAULT_SOURCE_EXT = .cpp

EXTRA_PROGRAMS = tests/tutorial tests/fib-example

tutorial: tests/tutorial
.PHONY: tutorial

fib-example: tests/fib-example
.PHONY: fib-example

pkgconfigdir=$(libdir)/pkgconfig
pkgconfig_DATA = libampp.pc
