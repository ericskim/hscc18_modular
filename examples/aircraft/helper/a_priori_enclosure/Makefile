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

TARGET = comp_enclosure 

all: $(TARGET)

%.o: %.cc
	$(CC) $(CXXFLAGS) $(SCOTSINC) $(VNODELPINC) -c -o $@ $<

comp_enclosure:  comp_enclosure.o
	$(CC) $(LDFLAGS) -o comp_enclosure  comp_enclosure.o $(VNODELPLIB) $(LIBS) 

clean:
	rm comp_enclosure.o comp_enclosure 
