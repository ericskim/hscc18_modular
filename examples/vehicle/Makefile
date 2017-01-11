#
# compiler
#
#CC        = g++
CC       	  = clang++
CXXFLAGS 		= -Wall -Wextra -std=c++11 -O3
CXXDEBUG		= -g -Wall -Wextra -std=c++11


#
# scots 
#
SCOTSROOT		= ../..
SCOTSINC		= -I$(SCOTSROOT)/sparse -I$(SCOTSROOT)/utils

SCOTSNEW		= -I$(SCOTSROOT)/src -I$(SCOTSROOT)/utils


.PHONY: vehicle vehicle_new

TARGET = vehicle vehicle_new

all: $(TARGET)

vehicle: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) vehicle.cc -o vehicle

vehicle_new: 
	$(CC) $(CXXFLAGS) $(SCOTSNEW) vehicle_new.cc -o vehicle_new


debug: 
	$(CC) $(CXXDEBUG) $(SCOTSNEW) vehicle_new.cc -o vehicle_new
	$(CC) $(CXXDEBUG) $(SCOTSINC) vehicle.cc -o vehicle


clean:
	rm vehicle 
