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


.PHONY: interconnection simulation interconnection_2 consensus consensus_simulation max_connection

TARGET = interconnection simulation interconnection_2 consensus consensus_simulation max_connection

all: $(TARGET)

interconnection: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) $(CUDD) interconnection.cc -o interconnection.o

interconnection_2: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) $(CUDD) interconnection_2.cc -o interconnection_2.o

consensus: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) $(CUDD) consensus.cc -o consensus.o

consensus_simulation: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) $(CUDD) consensus_simulation.cc -o consensus_simulation.o

simulation:
	$(CC) $(CXXFLAGS) $(SCOTSINC) $(CUDD) simulation.cc -o simulation.o

max_connection: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) $(CUDD) max_connection.cc -o max_connection.o

clean:
	rm interconnection simulation interconnection_2 consensus max_connection
