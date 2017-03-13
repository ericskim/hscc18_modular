#
# compiler
#
CXX       	  	= g++
#CXXFLAGS 		= -Wall -Wextra -std=c++11 -DSCOTS_BDD
CXXFLAGS 		= -Wall -Wextra -std=c++11 -O3 -DNDEBUG -DSCOTS_BDD
#
# scots 
#
SCOTSROOT		= ../..
SCOTSINC		= -I$(SCOTSROOT)/src  -I$(SCOTSROOT)/utils/
#
# cudd 
#
CUDDINC =  -I/usr/local/include 
CUDDLIB =  -L/usr/local/lib -lcudd


.PHONY: dcdc_bdd_recurrence 

TARGET = dcdc_bdd_recurrence 

all: $(TARGET)

dcdc_bdd_recurrence: 
	$(CXX) $(CXXFLAGS) $(CUDDINC) $(SCOTSINC)  dcdc_bdd_recurrence.cpp -o dcdc_bdd_recurrence $(CUDDLIB)

clean:
	rm dcdc_bdd_recurrence
