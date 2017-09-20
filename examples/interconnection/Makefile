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
CUDDPATH		=  /usr/local/Cellar/cudd/3.0.0
CUDDINC 		= -I$(CUDDPATH)/include
CUDD		    =  $(CUDDINC) -L$(CUDDPATH)/lib -lcudd


.PHONY: interconnection simulation interconnection_2

TARGET = interconnection simulation interconnection_2

all: $(TARGET)

interconnection: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) $(CUDD) interconnection.cc -o interconnection.o

interconnection_2: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) $(CUDD) interconnection_2.cc -o interconnection_2.o

simulation:
	$(CC) $(CXXFLAGS) $(SCOTSINC) $(CUDD) simulation.cc -o simulation.o

clean:
	rm interconnection simulation interconnection_2
