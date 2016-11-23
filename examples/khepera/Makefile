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
SCOTSINC		= -I$(SCOTSROOT)/src -I$(SCOTSROOT)/src/sparse


#
# cudd 
#


.PHONY: unicycle

TARGET = unicycle

all: $(TARGET)

unicycle: 
	$(CC) $(CXXFLAGS) $(SCOTSINC) unicycle.cc -o unicycle


clean:
	rm unicycle  *.o
