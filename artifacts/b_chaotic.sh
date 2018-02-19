#!/bin/sh
export BOOST_INSTALL=/home/users/p02119/software/install/boost
export AMPP_INSTALL=/home/users/p02119/software/install/ampp
export LIBCDS_INSTALL=/home/users/p02119/software/pbgl2/boost/external/pq/cds-2.1.0
export BOOST_ROOT=$BOOST_INSTALL

#fsanitize=address" \
rm -rf CMakeCache.txt CMakeFiles
make clean
cmake .. -DAM++_LIBRARY:FILEPATH=$AMPP_INSTALL/lib/libampp.a \
-DLIBCDS_INCLUDE_DIR:PATH=$LIBCDS_INSTALL \
-DLIBCDS_LIBRARY:FILEPATH=$LIBCDS_INSTALL/bin/gcc-amd64-linux-0/libcds.so \
-DCMAKE_EXE_LINKER_FLAGS:STRING="/home/users/p02119/software/pbgl2/boost/external/pq/cds-2.1.0/bin/gcc-amd64-linux-0/libcds.so -lnuma" \
-DLIBNBC_LIBRARY:FILEPATH=$LIBCDS_INSTALL/bin/gcc-amd64-linux-0/libcds.so \
-DAM++_INCLUDE_DIR:PATH=$AMPP_INSTALL/include \
-DLIBNBC_INCLUDE_DIR= \
-DCMAKE_CXX_FLAGS:STRING="-std=c++11 -ggdb -DDISABLE_SELF_SEND_CHECK -DAMPLUSPLUS_BUILTIN_ATOMICS -UDC_USE_PRIORITY -UDC_USE_EXPLICIT_POLLING -UDC_USE_HANDLERS_PENDING -UDC_USE_AGGRESSIVE_PRIORITY -UDC_USE_POLL_ONCE -fopenmp -fsanitize=address" \
-DBoost_INCLUDE_DIR:PATH=$BOOST_INSTALL/include \
-DBoost_LIBRARY_DIR:PATH=$BOOST_INSTALL/lib \
-DBoost_RANDOM_LIBRARY_DEBUG:FILEPATH=$BOOST_INSTALL/lib/libboost_random.a \
-DBoost_RANDOM_LIBRARY_RELEASE:FILEPATH=$BOOST_INSTALL/lib/libboost_random.a \
-DBoost_SYSTEM_LIBRARY_DEBUG:FILEPATH=$BOOST_INSTALL/lib/libboost_system.a \
-DBoost_SYSTEM_LIBRARY_RELEASE:FILEPATH=$BOOST_INSTALL/lib/libboost_system.a \
-DBoost_THREAD_LIBRARY_DEBUG:FILEPATH=$BOOST_INSTALL/lib/libboost_thread.a \
-DBoost_THREAD_LIBRARY_RELEASE:FILEPATH=$BOOST_INSTALL/lib/libboost_thread.a \
-DBoost_SYSTEM_LIBRARY:FILEPATH=$BOOST_INSTALL/lib/libboost_system.a \
-DBoost_THREAD_LIBRARY:FILEPATH=$BOOST_INSTALL/lib/libboost_thread.a \
-DCMAKE_CXX_COMPILER:FILEPATH=CC \
-DCMAKE_C_COMPILER:FILEPATH=cc \
-DMPI_LIBRARY:FILEPATH=/opt/cray/mpt/7.1.2/gni/mpich2-gnu/49/lib/libmpichcxx_gnu_49_mt.a  \
-DMPI_INCLUDE_PATH:PATH=/opt/cray/mpt/7.1.2/gni/mpich2-gnu/49/include


make performance_test VERBOSE=1
