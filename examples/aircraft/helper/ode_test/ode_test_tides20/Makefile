#
# compiler
#
#CC        = g++
C       	  = clang
CC       	  = clang++
CXXFLAGS 		= -g -Wall -Wextra -std=c++11
CXXFLAGS 		= -Wall -Wextra -std=c++11 -O3 -DNDEBUG
CFLAGS 		  = -Wall -Wextra -O3 


ROOT		= ../../..
#
# scots 
#
SCOTSINC		= -I$(ROOT)/sparse -I$(ROOT)/utils

#
# tide
#
TIDEPATH		= $(ROOT)/external/tides-2.0/libTIDES/


#
# cudd 
#


.PHONY: test

TARGET = test

all: $(TARGET)

%.o:%.c
	$(C)  $(CFLAGS) -I$(TIDEPATH) -c $< -o $@ 

%.o:%.cc
	$(CC) $(CXXFLAGS) $(SCOTSINC) -I$(TIDEPATH) -c $< -o $@ 

test: aircraftODE.o test.o
	$(CC) $(CXXFLAGS) -L$(TIDEPATH) -lm -lTIDES  aircraftODE.o test.o -o test


clean:
	rm test *.o
