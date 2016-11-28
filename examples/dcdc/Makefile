#
# compiler
#
#CC        = g++
CC       	  = g++
CXXFLAGS 		= -g -Wall -Wextra -std=c++14
CXXFLAGS 		= -Wall -Wextra -std=c++14 -O3 -DNDEBUG


#
# scots 
#
SCOTSROOT		= ../..
SCOTSINC		= -I$(SCOTSROOT)/sparse  -I$(SCOTSROOT)/utils/


#
# cudd 
#


.PHONY: dcdc

TARGET = dcdc

all: $(TARGET)

dcdc: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) dcdc.cc -o dcdc


clean:
	rm dcdc
