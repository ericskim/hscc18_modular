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
SCOTSINC		= -I$(SCOTSROOT)/sparse -I$(SCOTSROOT)/utils

SCOTSNEW		= -I$(SCOTSROOT)/src -I$(SCOTSROOT)/utils


#
# cudd 
#


.PHONY: aircraft aircraft_new

TARGET = aircraft aircraft_new

all: $(TARGET)

aircraft: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) aircraft.cc -o aircraft

aircraft_new: 
	$(CC) $(CXXFLAGS) $(SCOTSNEW) aircraft_new.cc -o aircraft_new


clean:
	rm aircraft aircraft_new
