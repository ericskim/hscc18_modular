#
# matlab root 
#
MATLABROOT	= /Applications/MATLAB_R2016b.app/

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
SCOTSINC		= -I$(SCOTSROOT)/src


all: sparse 

sparse: GridPoints mexStaticController


mexStaticController: mexStaticController.cc
	$(CC) COMPFLAGS='$(COMPFLAGS) $(CXXFLAGS)' -largeArrayDims $(SCOTSINC) mexStaticController.cc

GridPoints: GridPoints.cc
	$(CC) COMPFLAGS='$(COMPFLAGS) $(CXXFLAGS)' -largeArrayDims $(SCOTSINC) GridPoints.cc


