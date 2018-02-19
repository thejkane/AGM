BOOST_PATH=/Users/kane972/development/boost/install/
CXX=CC -std=c++11

INCLUDES=-I../include -I$(BOOST_PATH)/include -I.

ROOT_PATH=$(BOOST_PATH)/lib/

RMAT2=-DNEW_GRAPH500_SPEC

#LDFLAGS=-static -L$(BOOST_PATH)/lib/
LDFLAGS=-L$(BOOST_PATH)/lib/

CXXFILES = main.cpp

#CXXFLAGS = -O3 -o permute
CXXFLAGS =-Wall -ggdb -fno-omit-frame-pointer -dynamic -fsanitize=address -o gizmo

#LIBS = $(BOOST_PATH)/lib/libboost_system.a $(BOOST_PATH)/lib/libboost_thread.a
LIBS = 


all: 
	$(CXX) $(CXXFILES) $(LIBS) $(CXXFLAGS) $(INCLUDES)

clean: 
	rm -rf a.out gizmo *.o gizmo.dSYM

