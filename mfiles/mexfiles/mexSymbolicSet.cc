/*
 * mexSymbolicSet.cc
 *
 *  created: Jan 2017
 *   author: Matthias Rungger
 */

#include <iostream>
#include <vector>
/* mex */
#include "mex.h"
#include "ClassHandle.hh"
/* scots */
#include "scots.hh"



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  /* get pointer to the command */
  const char* command=mxArrayToString(prhs[0]);
  /* init */
  if(!strcmp(command,"init")) {
    /* get pointer to filename */
    const char* filename=mxArrayToString(prhs[1]);
    /* initialize cudd manager */ 
    Cudd* mgr = new Cudd;
    using input_type = std::array<double,1>;
    /* return a handle to the uniform grid instance */
    scots::SymbolicSet *set = new scots::SymbolicSet(*mgr,1,input_type{{.99}},input_type{{2.1}},input_type{{1}});
    BDD* bdd = new BDD;
    //scots::read_from_file(*set,*bdd,mgr,filename);
    //if(!scots::read_from_file(*set,*bdd,mgr,filename)) {
    //  plhs[0]=mxCreateDoubleMatrix(0,0,mxREAL);
    //  plhs[1]=mxCreateDoubleMatrix(0,0,mxREAL);
    //  return;
    //}
    plhs[0] = convertPtr2Mat<scots::SymbolicSet>(set);
    plhs[1] = convertPtr2Mat<BDD>(bdd);
    plhs[2] = convertPtr2Mat<Cudd>(mgr);
    return;
  }
  /* delete */
  if (!strcmp(command,"delete")) {
    /* delete the StaticController */
    scots::SymbolicSet *set= convertMat2Ptr<scots::SymbolicSet>(prhs[1]);
    set->print_info(1);
    destroyObject<scots::SymbolicSet>(prhs[1]);
    destroyObject<BDD>(prhs[2]);
    destroyObject<Cudd>(prhs[3]);
    return;
  }
  return;
}

