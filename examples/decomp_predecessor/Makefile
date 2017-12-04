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

 
.PHONY: interconnection simulation consensus consensus_simulation max_connection runningmax better_consensus iter_consensus monotone_consensus decomp_predecessor

TARGET = interconnection simulation consensus consensus_simulation max_connection runningmax better_consensus iter_consensus monotone_consensus decomp_predecessor

all: $(TARGET)

interconnection: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) $(CUDD) interconnection.cc -o interconnection.o

monotone_consensus:
	$(CC) $(CXXFLAGS) $(SCOTSINC) $(CUDD) monotone_consensus.cc -o monotone_consensus.o

consensus: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) $(CUDD) consensus.cc -o consensus.o

better_consensus:
	$(CC) $(CXXFLAGS) $(SCOTSINC) $(CUDD) better_consensus.cc -o better_consensus.o

iter_consensus:
	$(CC) $(CXXFLAGS) $(SCOTSINC) $(CUDD) iter_consensus.cc -o iter_consensus.o

consensus_simulation: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) $(CUDD) consensus_simulation.cc -o consensus_simulation.o

runningmax: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) $(CUDD) runningmax.cc -o runningmax.o

simulation:
	$(CC) $(CXXFLAGS) $(SCOTSINC) $(CUDD) simulation.cc -o simulation.o

max_connection: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) $(CUDD) max_connection.cc -o max_connection.o

decomp_predecessor: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) $(CUDD) decomp_predecessor.cc -o decomp_predecessor.o

clean:
	rm interconnection simulation interconnection_2 consensus max_connection decomp_predecessor
