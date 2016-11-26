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
SCOTSINC		= -I$(SCOTSROOT)/sparse\
							-I$(SCOTSROOT)/utils 


all: sparse 

sparse 		 	: mexUniformGrid\
							mexController


mexUniformGrid: mexUniformGrid.cc
	$(CC) COMPFLAGS='$(COMPFLAGS) $(CXXFLAGS)' -largeArrayDims $(SCOTSINC) mexUniformGrid.cc

mexController: mexController.cc
	$(CC) COMPFLAGS='$(COMPFLAGS) $(CXXFLAGS)' -largeArrayDims $(SCOTSINC) mexController.cc


