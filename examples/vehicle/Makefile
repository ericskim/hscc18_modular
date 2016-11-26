#
# compiler
#
#CC        = g++
CC       	  = clang++
CXXFLAGS 		= -Wall -Wextra -std=c++14 -O3  
CXXDEBUG		= -g -Wall -Wextra -std=c++14


#
# scots 
#
SCOTSROOT		= ../..
SCOTSINC		= -I$(SCOTSROOT)/sparse -I$(SCOTSROOT)/utils


.PHONY: vehicle

TARGET = vehicle

all: $(TARGET)

vehicle: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) vehicle.cc -o vehicle

debug: 
	$(CC) $(CXXDEBUG) $(SCOTSINC) vehicle.cc -o vehicle


clean:
	rm vehicle  *.o
