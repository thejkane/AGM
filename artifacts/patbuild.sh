#!/bin/sh

module load perftools-base
module load perftools

export BINARY=tc_family

export BOOST_INSTALL=/home/users/p02119/development/install/boost/
export AMPP_INSTALL=/home/users/p02119/development/install/amppvar/
export LIBCDS_INSTALL=/home/users/p02119/development/cds-2.1.0/
export BOOST_ROOT=$BOOST_INSTALL

export CRAYPE_LINK_TYPE=dynamic
#====================================================================================================================================#
binfile="./bin/$BINARY"
if [ -f "$binfile" ]
then
	rm $binfile
fi

rm -rf CMakeCache.txt CMakeFiles
make clean
cmake .. -D-DCMAKE_BUILD_TYPE=Release -DAM++_LIBRARY:FILEPATH=$AMPP_INSTALL/lib/libampp.a \
-DLIBCDS_INCLUDE_DIR:PATH=$LIBCDS_INSTALL \
-DLIBCDS_LIBRARY:FILEPATH=/home/users/p02119/development/cds-2.1.0/bin/cc-amd64-linux-64/libcds.so \
-DAM++_INCLUDE_DIR:PATH=$AMPP_INSTALL/include \
-DLIBNBC_INCLUDE_DIR= \
-DCMAKE_CXX_FLAGS:STRING="-g -mavx -march=corei7-avx -mtune=corei7-avx -std=c++11 -O3 -Ofast -ffast-math -funroll-loops -UAMPLUSPLUS_ENABLE_PERFORMANCE_COUNTERS -DAMPLUSPLUS_BUILTIN_ATOMICS -UDC_USE_PRIORITY -UDC_USE_EXPLICIT_POLLING -UDC_USE_HANDLERS_PENDING -UDC_USE_AGGRESSIVE_PRIORITY -UDC_USE_POLL_ONCE -DCYCLIC_DIST -DCC_PRIORITY -DTC_STATS -DTC_NO_ORDERING" \
-DCMAKE_EXE_LINKER_FLAGS:STRING="-L/home/users/p02119/development/cds-2.1.0/bin/cc-amd64-linux-64/ -lcds -lnuma" \
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
-DCMAKE_C_COMPILER:FILEPATH=cc

make tc_family VERBOSE=1

outfile="./bin/$BINARY"
if [ -f "$outfile" ]
then
	/bin/mv $outfile $binfile
else
	echo "$outfile not found."
	exit -1
fi

#====================================================================================================================================#






pat_build -w -g mpi $binfile
#pat_build -u -Dtrace-text-size=800 $binfile
#pat_build -S $binfile
mv *+pat ./bin/
