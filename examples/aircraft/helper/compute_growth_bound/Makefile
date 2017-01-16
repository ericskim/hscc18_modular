#
# compiler
#
CC       	  = clang++
CXXFLAGS 		= -O3 -Wextra -Wall -std=c++14 -DPROFIL_VNODE -DMAXORDER=50

LDFLAGS 		= 

#
# scots 
#
SCOTSROOT		= ../../../..
SCOTSINC		= -I$(SCOTSROOT)/src -I$(SCOTSROOT)/utils

#
# vnodelp/fadbad/filib 
#
VNODELPINC  = -I$(SCOTSROOT)/external/Profil-2.0.8/include \
              -I$(SCOTSROOT)/external/Profil-2.0.8/include/BIAS/ \
              -I$(SCOTSROOT)/external/vnodelp/FADBAD++ \
              -I$(SCOTSROOT)/external/vnodelp/include

VNODELPLIB = -L$(SCOTSROOT)/external/vnodelp/lib \
             -L$(SCOTSROOT)/external/Profil-2.0.8/lib \
             -L$(SCOTSROOT)/external/usr/local/lib 

LIBS = -lvnode -lProfil -lBias -llr

TARGET = comp_growth_bound comp_growth_bound_per_input

all: $(TARGET)

%.o: %.cc
	$(CC) $(CXXFLAGS) $(SCOTSINC) $(VNODELPINC) -c -o $@ $<

comp_growth_bound:  comp_growth_bound.o
	$(CC) $(LDFLAGS) -o comp_growth_bound  comp_growth_bound.o $(VNODELPLIB) $(LIBS) 

comp_growth_bound_per_input:  comp_growth_bound_per_input.o
	$(CC) $(LDFLAGS) -o comp_growth_bound_per_input  comp_growth_bound_per_input.o $(VNODELPLIB) $(LIBS) 


clean:
	rm comp_growth_bound.o comp_growth_bound comp_growth_bound_per_input.o comp_growth_bound_per_input
