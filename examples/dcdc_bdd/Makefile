#
# compiler
#
#CC        = g++
CC       	  = clang++
CXXFLAGS 		= -g -Wall -Wextra -std=c++11 -DSCOTS_BDD
CXXFLAGS 		= -Wall -Wextra -std=c++11 -O3 -DNDEBUG -DSCOTS_BDD
#
# scots 
#
SCOTSROOT		= ../..
SCOTSINC		= -I$(SCOTSROOT)/src  -I$(SCOTSROOT)/utils/
#
# cudd 
#
CUDD		    =  -I/opt/local/include -L/opt/local/lib -lcudd


.PHONY: dcdc simulate

TARGET = dcdc simulate

all: $(TARGET)

dcdc: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) $(CUDD) dcdc.cc -o dcdc

simulate: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) $(CUDD) simulate.cc -o simulate

clean:
	rm dcdc simulate
