#
# compiler
#
#CC        = g++
C       	  = clang
CC       	  = clang++
CXXFLAGS 		= -g -Wall -Wextra -std=c++11
CXXFLAGS 		= -Wall -Wextra -std=c++11 -O3 -DNDEBUG
CFLAGS 		= -Wall -Wextra -O3 


#
# scots 
#
SCOTSROOT		= ../..
SCOTSINC		= -I$(SCOTSROOT)/src -I/opt/local/include

#
# tide
#
TIDEPATH		= ../../../tides-2.0/libTIDES/


#
# cudd 
#


.PHONY: vehicle

TARGET = vehicle

all: $(TARGET)

%.o:%.c
	$(C)  $(CFLAGS) -I$(TIDEPATH) -c $< -o $@ 

%.o:%.cc
	$(CC) $(CXXFLAGS) $(SCOTSINC) -I$(TIDEPATH) -c $< -o $@ 

vehicle: vehicleODE.o vehicle.o
	$(CC) $(CXXFLAGS) -L$(TIDEPATH) -lm -lTIDES  vehicleODE.o vehicle.o -o vehicle


clean:
	rm vehicle  *.o
