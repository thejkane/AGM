#!/bin/sh
export BOOST_INSTALL=/home/users/p02119/development/install/boost/
export AMPP_INSTALL=/home/users/p02119/development/install/amppvar/
export BOOST_ROOT=$BOOST_INSTALL

cd ../
make install
cd tests

#export INCS='-I$BOOST_INSTALL/include -I$AMPP_INSTALL/include'
#export LIBS='$BOOST_INSTALL/lib/libboost_random.a $BOOST_INSTALL/lib/libboost_random.a $BOOST_INSTALL/lib/libboost_system.a $BOOST_INSTALL/lib/libboost_thread.a'

#CC -dynamic -fsanitize=address -fsanitize=signed-integer-overflow -std=c++11 testbuffer.cpp -I$BOOST_INSTALL/include -I$AMPP_INSTALL/include -I../ $BOOST_INSTALL/lib/libboost_random.a $BOOST_INSTALL/lib/libboost_system.a ../.libs/libampp.a $BOOST_INSTALL/lib/libboost_thread.a -o varbuffer -v
CC -dynamic -fsanitize=address -fsanitize=signed-integer-overflow -std=c++11 testsizebuffer.cpp -I$BOOST_INSTALL/include -I$AMPP_INSTALL/include -I../ $BOOST_INSTALL/lib/libboost_random.a $BOOST_INSTALL/lib/libboost_system.a ../.libs/libampp.a $BOOST_INSTALL/lib/libboost_thread.a -o varbuffer -v -UAMPLUSPLUS_USE_THREAD_SERIALIZED -UNDEBUG
