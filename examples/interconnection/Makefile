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


.PHONY: consensus runningmax

TARGET = consensus runningmax

all: $(TARGET)

consensus:
	$(CC) $(CXXFLAGS) $(SCOTSINC) $(CUDD) consensus.cc -o consensus.o

runningmax: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) $(CUDD) runningmax.cc -o runningmax.o

clean:
	rm consensus runningmax
