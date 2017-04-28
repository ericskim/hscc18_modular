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


.PHONY: sparse_test sparse_vehicle

TARGET = sparse_test sparse_vehicle
#simulate

all: $(TARGET)

sparse_test: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) $(CUDD) sparse_test.cc -o sparse_test

sparse_vehicle: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) $(CUDD) sparse_vehicle.cc -o sparse_vehicle

clean:
	rm sparse_test sparse_vehicle
