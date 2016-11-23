#
# compiler
#
#CC        = g++
CC       	  = g++
CXXFLAGS 		= -g -Wall -Wextra -std=c++11
CXXFLAGS 		= -g -Wall -Wextra -std=c++11
CXXFLAGS 		= -Wall -Wextra -std=c++11 -O3 -DNDEBUG


#
# scots 
#
SCOTSROOT		= ../..
SCOTSINC		= -I$(SCOTSROOT)/src/sparse  -I$(SCOTSROOT)/src/


#
# cudd 
#


.PHONY: dcdc

TARGET = dcdc

all: $(TARGET)

dcdc: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) dcdc.cc -o dcdc


clean:
	rm dcdc  *.o
