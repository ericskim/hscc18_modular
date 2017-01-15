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
SCOTSINC		= -I$(SCOTSROOT)/src  -I$(SCOTSROOT)/utils/

.PHONY: dcdc simulate

TARGET = dcdc simulate

all: $(TARGET)

dcdc: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) dcdc.cc -o dcdc

simulate: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) simulate.cc -o simulate

clean:
	rm dcdc simulate
