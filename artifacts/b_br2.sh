#!/bin/sh
export BOOST_INSTALL=/N/u/thejkane/BigRed2/crest/development/install/boost
export AMPP_INSTALL=/N/u/thejkane/BigRed2/crest/development/install/ampp
export LIBCDS_INSTALL=/N/dc2/scratch/thejkane/cds-1.6.0/
export BOOST_ROOT=$BOOST_INSTALL

rm -rf CMakeCache.txt CMakeFiles
make clean
cmake .. -DAM++_LIBRARY:FILEPATH=$AMPP_INSTALL/lib/libampp.a \
-DLIBCDS_INCLUDE_DIR:PATH=$LIBCDS_INSTALL \
-DLIBCDS_LIBRARY:FILEPATH=/N/dc2/scratch/thejkane/cds-1.6.0/bin/gcc-amd64-linux-0/libcds.so \
-DAM++_INCLUDE_DIR:PATH=$AMPP_INSTALL/include \
-DLIBNBC_INCLUDE_DIR= \
-DCMAKE_CXX_FLAGS:STRING="-std=c++11 -ggdb -DAMPLUSPLUS_BUILTIN_ATOMICS -UDC_USE_PRIORITY -UDC_USE_EXPLICIT_POLLING -UDC_USE_HANDLERS_PENDING -UDC_USE_AGGRESSIVE_PRIORITY -UDC_USE_POLL_ONCE" \
-DCMAKE_EXE_LINKER_FLAGS:STRING="-dynamic -L/N/dc2/scratch/thejkane/cds-1.6.0/bin/gcc-amd64-linux-0/ -lcds -lnuma" \
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
-DMPI_LIBRARY:FILEPATH=/opt/cray/mpt/7.2.6/gni/mpich-gnu/51/lib/libmpichcxx_gnu_51_mt.a  \
-DMPI_INCLUDE_PATH:PATH=/opt/cray/mpt/7.2.6/gni/mpich-gnu/51/include/


make performance_test VERBOSE=1
