#
# compiler
#
#CC        = g++
CC       	  = clang++
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


.PHONY: vehicle simulate

TARGET = vehicle simulate

all: $(TARGET)

vehicle: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) $(CUDD) vehicle.cc -o vehicle

simulate: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) $(CUDD) simulate.cc -o simulate

clean:
	rm vehicle simulate
