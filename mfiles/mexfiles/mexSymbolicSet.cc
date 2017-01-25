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
    /* initialize cudd manager and load BDD and SymbolicSet from file */ 
    Cudd* mgr = new Cudd;
    BDD* bdd = new BDD;
    scots::SymbolicSet *set = new scots::SymbolicSet;
    if(!scots::read_from_file(*set,*bdd,*mgr,filename)) {
      plhs[0]=mxCreateDoubleMatrix(0,0,mxREAL);
      plhs[1]=mxCreateDoubleMatrix(0,0,mxREAL);
      return;
    }
    /* store ptr into return value */
    plhs[0]=mxCreateNumericMatrix(3,1, mxUINT64_CLASS, mxREAL);
    uint64_t *p=(uint64_t*)mxGetData(plhs[0]);
    p[0]=convertPtr2MatUINT64<scots::SymbolicSet>(set);
    p[1]=convertPtr2MatUINT64<BDD>(bdd);
    p[2]=convertPtr2MatUINT64<Cudd>(mgr);
    /* dimension */
    plhs[1]=mxCreateDoubleMatrix(1,1, mxREAL);
    double *d=mxGetPr(plhs[1]);
    d[0]=set->get_dim();
    return;
  }

  if(mxGetClassID(prhs[1])!=mxUINT64_CLASS)
    mexErrMsgTxt("Second argument to mexSymbolicSet must be a real uint64 array.");
  uint64_t* p=(uint64_t*) mxGetData(prhs[1]);

  /* delete */
  if (!strcmp(command,"delete")) {
    /* delete the StaticController */
    destroyObject<scots::SymbolicSet>(p[0]);
    destroyObject<BDD>(p[1]);
    destroyObject<Cudd>(p[2]);
    return;
  }

  /* extract pointer information */
  scots::SymbolicSet *set = convertMat2Ptr<scots::SymbolicSet>(p[0]);
  BDD* bdd = convertMat2Ptr<BDD>(p[1]);
  Cudd* mgr = convertMat2Ptr<Cudd>(p[2]);

  /* return grid points */ 
  if(!strcmp(command,"gridpoints")) {
    /* get projection indices */
    double* indices=mxGetPr(prhs[2]);
    const mwSize *dims=mxGetDimensions(prhs[2]);
    std::vector<int> project_dim {}; 
    for(int i=0; i<dims[0]; i++) 
      project_dim.push_back(indices[i]);

    auto grid_points = set->projection(*mgr,*bdd,project_dim);
    size_t no_grid_points = grid_points.size()/dims[0];
    /* convert abstract cell IDs to the center of cells */
    plhs[0]=mxCreateDoubleMatrix(no_grid_points,dims[0],mxREAL);
    double *ptr=mxGetPr(plhs[0]);
    for(size_t i=0; i<no_grid_points; i++) {
      for(int j=0; j<dims[0]; j++) {
        ptr[j*no_grid_points+i]=grid_points[i*dims[0]+j];
      }
    }
    return;
  } 

  /* return grid points restricted to x in rhs[1] */ 
  if(!strcmp(command,"restriction")) {
    /* get restriction state */
    double* xp=mxGetPr(prhs[2]);
    double* ip=mxGetPr(prhs[3]);
    const mwSize *dims=mxGetDimensions(prhs[2]);
    std::vector<double> xv {}; 
    std::vector<int> iv {}; 
    for(int i=0; i<dims[0]; i++) {
      xv.push_back(xp[i]);
      iv.push_back(ip[i]);
    }
    size_t codim = set->get_dim()-dims[0];
    auto grid_points = set->restriction(*mgr,*bdd,xv,iv);
    if(!grid_points.size() || !codim) {
      plhs[0]=mxCreateDoubleMatrix(0,0,mxREAL);
      return;
    }
    size_t no_grid_points = grid_points.size()/codim;
    /* convert abstract cell IDs to the center of cells */
    plhs[0]=mxCreateDoubleMatrix(no_grid_points,codim,mxREAL);
    double *ptr=mxGetPr(plhs[0]);
    for(size_t i=0; i<no_grid_points; i++) {
      for(int j=0; j<codim; j++) {
        ptr[j*no_grid_points+i]=grid_points[i*codim+j];
      }
    }
    return;
  } 
  return;
}

