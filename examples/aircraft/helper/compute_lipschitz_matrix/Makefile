#
# compiler
#
CC       	  = clang++
CXXFLAGS 		= -O3 -Wextra -Wall -std=c++11 -DPROFIL_VNODE -DMAXORDER=50

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

TARGET = comp_lipschitz_matrix

all: $(TARGET)

%.o: %.cc
	$(CC) $(CXXFLAGS) $(SCOTSINC) $(VNODELPINC) -c -o $@ $<

comp_lipschitz_matrix:  comp_lipschitz_matrix.o
	$(CC) $(LDFLAGS) -o comp_lipschitz_matrix  comp_lipschitz_matrix.o $(VNODELPLIB) $(LIBS) 

clean:
	rm comp_lipschitz_matrix.o comp_lipschitz_matrix 
