#
# matlab root 
#
MATLABROOT	= /Applications/MATLAB_R2016a.app/

#
# compiler
#
#CC        = g++
CC       	  = $(MATLABROOT)/bin/mex
CXXFLAGS 		= -g -Wall -Wextra -std=c++11 -O3 

#
# scots
#
SCOTSROOT   = ../..
SCOTSINC		= -I$(SCOTSROOT)/src/\
							-I$(SCOTSROOT)/src/sparse \
							-I$(SCOTSROOT)/src/bdd\
							-I$(SCOTSROOT)/src/linkedlist

#
# cudd 
#
CUDDPATH		=  /opt/local/
CUDDINC 		= -I$(CUDDPATH)/include
CUDDLIBS		= -lcudd 
CUDDLPATH   = -L$(CUDDPATH)/lib


all: sparse bdd

sparse 		 	: mexUniformGrid\
							mexController

bdd 			: mexSymbolicSetReadFromFile


mexSymbolicSetReadFromFile: mexSymbolicSetReadFromFile.cc
	$(CC) COMPFLAGS='$(COMPFLAGS) $(CXXFLAGS)' -largeArrayDims $(SCOTSINC) $(CUDDINC) $(CUDDLPATH) $(CUDDLIBS) mexSymbolicSetReadFromFile.cc

mexUniformGrid: mexUniformGrid.cc
	$(CC) COMPFLAGS='$(COMPFLAGS) $(CXXFLAGS)' -largeArrayDims $(SCOTSINC) mexUniformGrid.cc

mexController: mexController.cc
	$(CC) COMPFLAGS='$(COMPFLAGS) $(CXXFLAGS)' -largeArrayDims $(SCOTSINC) mexController.cc


