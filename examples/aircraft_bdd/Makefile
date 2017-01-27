#
# compiler
#
#CC        = g++
CC       	  = clang++
CXXFLAGS 		= -g -Wall -Wextra -std=c++11
CXXFLAGS 		= -Wall -Wextra -std=c++11 -O3 -DNDEBUG -DSCOTS_BDD
#
# scots 
#
SCOTSROOT		= ../..
SCOTSINC		= -I$(SCOTSROOT)/src -I$(SCOTSROOT)/utils
#
# cudd 
#
CUDD		    =  -I/opt/local/include -L/opt/local/lib -lcudd

.PHONY: aircraft simulate

TARGET = aircraft simulate

all: $(TARGET)

aircraft: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) $(CUDD) aircraft.cc -o aircraft

simulate: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) $(CUDD) simulate.cc -o simulate


clean:
	rm aircraft simulate
