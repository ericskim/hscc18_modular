#
# matlab root 
#
MATLABROOT	= /Applications/MATLAB_R2016b.app/
#
# compiler
#
CC       	  = $(MATLABROOT)/bin/mex
CXXFLAGS 		= -g -Wall -Wextra -std=c++11 -O3 
#
# scots
#
SCOTSROOT   = ../..
SCOTSINC		= -I$(SCOTSROOT)/src
#
# cudd 
#
CUDDPATH		=  /opt/local/
CUDDINC 		= -I$(CUDDPATH)/include
CUDDLIBS		= -lcudd 
CUDDLPATH   = -L$(CUDDPATH)/lib



all: sparse bdd

sparse: GridPoints mexStaticController

bdd: mexSymbolicSet


mexStaticController: mexStaticController.cc
	$(CC) COMPFLAGS='$(COMPFLAGS) $(CXXFLAGS)' -largeArrayDims $(SCOTSINC) mexStaticController.cc

GridPoints: GridPoints.cc
	$(CC) COMPFLAGS='$(COMPFLAGS) $(CXXFLAGS)' -largeArrayDims $(SCOTSINC) GridPoints.cc

mexSymbolicSet: mexSymbolicSet.cc
	$(CC) COMPFLAGS='$(COMPFLAGS) $(CXXFLAGS)' -DSCOTS_BDD -largeArrayDims $(SCOTSINC) $(CUDDINC) $(CUDDLPATH) $(CUDDLIBS) mexSymbolicSet.cc

