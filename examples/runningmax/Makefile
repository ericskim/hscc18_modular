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

#
# example flags
#
EXFLAGS 		= -D NUMSYS=5

.PHONY: runningmax

TARGET = runningmax

all: $(TARGET)

runningmax: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) $(CUDD) $(EXFLAGS) runningmax.cc -o runningmax.o

clean:
	rm runningmax
