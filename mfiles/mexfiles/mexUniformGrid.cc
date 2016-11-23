/*
 * mexUniformGrid.cc
 *
 *  created on: 11.01.2016
 *      author: rungger
 */


#include <iostream>
#include <math.h>
#include <cstring>
#include <stdexcept> 


#include "mex.h"
#include "ClassHandle.hh"
#include "UniformGrid.hh"
#include "IO.hh"

typedef scots::UniformGrid<std::vector<double> > Grid;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  /* there must be at least and it must be a string */
  if (nrhs < 1 || mxIsChar(prhs[0]) != 1)
    mexErrMsgIdAndTxt("MATLAB:mexUniformGrid","mexUniformGrid requires at least one input argument.");
  /* get pointer to the command */
  const char* command=mxArrayToString(prhs[0]);
  /* init */
  if(!strcmp(command,"init")) {
    /* there must a second parameter containing the file name */
    if (nrhs < 2 || mxIsChar(prhs[1]) != 1)
      mexErrMsgIdAndTxt("MATLAB:mexUniformGrid","'init' requires the filename where the uniform grid is saved as second parameter.");
    /* there must be one return parameter */
    if (nlhs != 1)
        mexErrMsgTxt("MATLAB:mexUniformGrid.hh: There must be one return parameter.");
    const char* filename=mxArrayToString(prhs[1]);
    /* return a handle to the uniform grid instance */
    Grid* gs=new Grid;
    scots::IO::readFromFile(gs,filename);
    plhs[0] = convertPtr2Mat<Grid>(gs);
    return;
  }
  /* delete */
  if (!strcmp(command,"delete")) {
    /* free the uniform grid object */
    destroyObject<Grid>(prhs[1]);
    return;
  }

  if(nrhs<2)
    mexErrMsgTxt("MATLAB:mexUniformGrid.hh: Second input should be the object handle.");
  /* get pointer to UniformGrid object */ 
  Grid *grid = convertMat2Ptr<Grid>(prhs[1]);
  /* get state space dimension */
  size_t dim=grid->getDimension();

  /* set: return grid points in abstract set  */
  if (!strcmp(command,"set")) {
    if(nrhs>2 && mxIsNumeric(prhs[2])) {
      mwSize n=mxGetN(prhs[2]);
      double* ptr=mxGetPr(prhs[2]);
      std::vector<size_t> projectdim;
      for(size_t i=0; i<n; i++) {
        projectdim.push_back((size_t)ptr[i]);
      }
      /* check sanity of project dimensions */
      if(n> dim)
        mexErrMsgIdAndTxt( "MATLAB:mexUniformGrid", "Number of project dimensions cannot exceed the number of dimensions.");
      for(size_t i=0; i<n; i++) {
        if(projectdim[i]>dim)
          mexErrMsgIdAndTxt( "MATLAB:mexUniformGrid", "Cannot project UniformGrid onto given dimensions.");
      }
      /* create matrix to store grid points in the projected abstract set */
      std::vector<size_t> setProj=grid->projectSet(projectdim);
      size_t size=setProj.size();
      if(size) {
        plhs[0]=mxCreateDoubleMatrix((mwSize)size,(mwSize)n,mxREAL);
        double *points=mxGetPr(plhs[0]);
        for(size_t v=0;v<size;v++) {
          std::vector<double> x(dim);
          grid->itox(setProj[v],x);
          for(size_t w=0;w<n;w++)
            points[w*size+v]=x[projectdim[w]-1];
        }
      } else
        plhs[0]=mxCreateDoubleScalar(-1);
      return;
    }
    else if(nrhs==2) {
      /* create matrix to store grid points in the abstract set */
      size_t size=grid->getSetSize();
      if(size) {
        plhs[0]=mxCreateDoubleMatrix((mwSize)(dim*size),1,mxREAL);
        double *points=mxGetPr(plhs[0]);
        grid->setIdxToDoubleMatrix(points);
        mwSize ndim[2]={size, dim};
        mxSetDimensions(plhs[0],ndim,2);
      } else
        plhs[0]=mxCreateDoubleScalar(-1);
      return;
    }
  }
    
  /* dim */
  if (!strcmp(command,"dim")) {
    plhs[0]=mxCreateDoubleScalar((double)dim);
  }
  /* eta */
  if (!strcmp(command,"eta")) {
    std::vector<double> eta(dim);
    eta = grid->getEta();
    plhs[0]=mxCreateDoubleMatrix((mwSize)dim,1,mxREAL);
    std::copy(eta.begin(), eta.end(), mxGetPr(plhs[0]));
  }
  /* first */
  if (!strcmp(command,"first")) {
    std::vector<double> first(dim);
    first = grid->getFirstGridPoint();
    plhs[0]=mxCreateDoubleMatrix((mwSize)dim,1,mxREAL);
    std::copy(first.begin(), first.end(), mxGetPr(plhs[0]));
  }
  /* z - measurement error bound */
  if (!strcmp(command,"zzz")) {
    std::vector<double> zzz(dim);
    zzz = grid->getZ();
    plhs[0]=mxCreateDoubleMatrix((mwSize)dim,1,mxREAL);
    std::copy(zzz.begin(), zzz.end(), mxGetPr(plhs[0]));
  }
  /* nofgp */
  if (!strcmp(command,"nofgp")) {
    std::vector<size_t> nofgp(dim);
    nofgp = grid->getNofGridPoints();
    plhs[0]=mxCreateDoubleMatrix((mwSize)dim,1,mxREAL);
    std::copy(nofgp.begin(), nofgp.end(), mxGetPr(plhs[0]));
  }
  return;
}

