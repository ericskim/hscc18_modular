#
# compiler
#
#CC        = g++
CC       	  = clang++
CXXFLAGS 		= -g -Wall -Wextra -std=c++11
CXXFLAGS 		= -Wall -Wextra -std=c++11 -O3 -DNDEBUG


#
# scots 
#
SCOTSROOT		= ../..
SCOTSINC		= -I$(SCOTSROOT)/src -I$(SCOTSROOT)/utils


.PHONY: aircraft simulate

TARGET = aircraft simulate

all: $(TARGET)

aircraft: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) aircraft.cc -o aircraft

simulate: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) simulate.cc -o simulate


clean:
	rm aircraft simulate
