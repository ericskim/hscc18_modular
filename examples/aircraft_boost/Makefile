#
# compiler
#
#CC        = g++
CC       	  = clang++
CXXFLAGS 		= -g -Wall -Wextra -std=c++11
CXXFLAGS 		= -Wall -Wextra -std=c++11 -O3 
#
# scots 
#
SCOTSROOT		= ../..
SCOTSINC		= -I$(SCOTSROOT)/src -I$(SCOTSROOT)/utils
#
# boost 
#
BOOSTINC		= -I/opt/local/include


.PHONY: aircraft simulate

TARGET = aircraft simulate

all: $(TARGET)

aircraft: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) $(BOOSTINC) aircraft.cc -o aircraft

simulate: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) $(BOOSTINC) simulate.cc -o simulate


clean:
	rm aircraft simulate
