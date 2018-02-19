BOOST_PATH=/home/users/p02119/development/install/boost/
AMPP_INSTALL=/home/users/p02119/development/install/amppvar/
CXX=CC -std=c++11

INCLUDES=-I$(AMPP_INSTALL)/include -I$(BOOST_PATH)/include 


#LDFLAGS=-static -L$(BOOST_PATH)/lib/
LDFLAGS=-L$(BOOST_PATH)/lib/

CXXFILES = latency_calculator.cpp

#CXXFLAGS = -DCYCLE_TYPE -O3 -o calc 
CXXFLAGS = -O3 -o calc 
#CXXFLAGS =-Wall -ggdb -fno-omit-frame-pointer -dynamic -fsanitize=address -o calc

#LIBS = $(BOOST_PATH)/lib/libboost_system.a $(BOOST_PATH)/lib/libboost_thread.a
LIBS = $(AMPP_INSTALL)/lib/libampp.a $(BOOST_PATH)/lib/libboost_system.a \
	$(BOOST_PATH)/lib/libboost_thread.a


all: 
	$(CXX) $(CXXFILES) $(LIBS) $(CXXFLAGS) $(INCLUDES)

clean: 
	rm -rf a.out calc

