#!/bin/sh
export BOOST_INSTALL=/Users/kane972/development/boost/install_1_55/
export AMPP_INSTALL=/Users/kane972/development/repos/install/ampp/
export LIBCDS_INSTALL=/Users/kane972/development/cds/cds-2.1.0/
export C_COMPILER=mpicc
export CXX_COMPILER=mpicxx
export BOOST_ROOT=$BOOST_INSTALL
export CRAYPE_LINK_TYPE=dynamic

rm -rf CMakeCache.txt CMakeFiles
make clean
cmake .. -D-DCMAKE_BUILD_TYPE=Release -DAM++_LIBRARY:FILEPATH=$AMPP_INSTALL/lib/libampp.a \
-DLIBCDS_INCLUDE_DIR:PATH=$LIBCDS_INSTALL \
-DLIBCDS_LIBRARY:FILEPATH=/Users/kane972/development/cds/cds-2.1.0/bin/gcc-amd64-darwin-64/libcds.dylib \
-DAM++_INCLUDE_DIR:PATH=$AMPP_INSTALL/include \
-DLIBNBC_INCLUDE_DIR= \
-DCMAKE_CXX_FLAGS:STRING="-std=c++11 -O3 -UAMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS -DAMPLUSPLUS_BUILTIN_ATOMICS -UDC_USE_PRIORITY -UDC_USE_EXPLICIT_POLLING -UDC_USE_HANDLERS_PENDING -UDC_USE_AGGRESSIVE_PRIORITY -UDC_USE_POLL_ONCE -DAGM_STATS -DDISABLE_SELF_SEND_CHECK -DNO_COALESCING -DPBGL2_PRINT_WORK_STATS" \
-DCMAKE_EXE_LINKER_FLAGS:STRING="-L/Users/kane972/development/cds/cds-2.1.0/bin/gcc-amd64-darwin-64/ -lcds" \
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
-DCMAKE_CXX_COMPILER:FILEPATH=$CXX_COMPILER \
-DCMAKE_C_COMPILER:FILEPATH=$C_COMPILER

make $1 VERBOSE=1

