 /*
 * mexReachController.cc
 *
 *  created on: 11.01.2016
 *      author: rungger
 */

#include <iostream>
#include <math.h>
#include <cstring>
#include <stdexcept> 

#include "mex.h"
#include "IO.hh"
#include "ClassHandle.hh"
#include "UniformGrid.hh"
#include "Game.hh"

typedef scots::Game Game;
typedef scots::UniformGrid<std::vector<double> > Grid;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  /* there must be at least and it must be a string */
  if (nrhs < 1 || mxIsChar(prhs[0]) != 1)
    mexErrMsgIdAndTxt("MATLAB:mexReachController","mexReachController requires at least one input argument.");
  /* get pointer to the command */
  const char* command=mxArrayToString(prhs[0]);
  /* init */
  if(!strcmp(command,"init")) {
    /* there must a second parameter containing the file name */
    if (nrhs < 2 || mxIsChar(prhs[1]) != 1)
      mexErrMsgIdAndTxt("MATLAB:mexReachController","'init' requires the filename where the uniform grid is saved as second parameter.");
    /* there must be three return parameter */
    if (nlhs != 3)
        mexErrMsgTxt("MATLAB:mexReachController.hh: There must be three return parameter (one for each object handle).");
    const char* filename=mxArrayToString(prhs[1]);
    /* return a handle to the uniform grid instance */
    Game* gg=new Game;
    Grid* stateGrid=new Grid;
    Grid* inputGrid=new Grid;
    scots::IO::readControllerFromFile(gg,filename);
    scots::IO::readFromFile(stateGrid,filename,"state");
    scots::IO::readFromFile(inputGrid,filename,"input");
    plhs[0] = convertPtr2Mat<Game>(gg);
    /* return a handle to the state space uniform grid instance */
    plhs[1] = convertPtr2Mat<Grid>(stateGrid);
    /* return a handle to the input space uniform grid instance */
    plhs[2] = convertPtr2Mat<Grid>(inputGrid);

    return;
  }
  /* delete */
  if (!strcmp(command,"delete")) {
    if (nrhs < 4)
      mexErrMsgIdAndTxt("MATLAB:mexReachController","'delete' requires four input parameters.");
    /* delete the UniformGrid and ReachabilityGame object */
    destroyObject<Game>(prhs[1]);
    destroyObject<Grid>(prhs[2]);
    destroyObject<Grid>(prhs[3]);
    return;
  }

  if(nrhs<4)
    mexErrMsgTxt("MATLAB:mexUniformGrid.hh: Four inputs expected.");
  /* get pointer to UniformGrid object */ 
  Game *game = convertMat2Ptr<Game>(prhs[1]);
  Grid *ss = convertMat2Ptr<Grid>(prhs[2]);
  Grid *is = convertMat2Ptr<Grid>(prhs[3]);

  /* state space dimension */
  size_t sdim=ss->getDimension();
  size_t idim=is->getDimension();

  /* input: return control input associated with x */
  if (!strcmp(command,"input")) {
    double *p=mxGetPr(prhs[4]);
    size_t xi;
    std::vector<size_t> ui;
    std::vector<size_t> uiProj;
    std::vector<double> x(sdim);
    std::vector<double> u(idim);
    /* copy state to std vector */
    for(size_t i=0; i<sdim; i++)
      x[i]=p[i];
    /* get index associated with  x */
    ss->xtoi(xi,x);
    /* get input idx associated with state index xi */
    ui = game->getLabel(xi);
    if(ui.size()==0) {
      plhs[0]=mxCreateDoubleMatrix(0,0,mxREAL);
      return;
    }


    /* get input associated with index ui */
    /* create matrix to store input */
    if(nrhs>5 && mxIsNumeric(prhs[5])) {
      mwSize n=mxGetN(prhs[5]);
      double *ptr=mxGetPr(prhs[5]);
      std::vector<size_t> projectdim;
      for(size_t i=0; i<n; i++) {
        projectdim.push_back((size_t)ptr[i]);
      }
      if(n> idim)
        mexErrMsgIdAndTxt( "MATLAB:mexController", "Number of project dimensions cannot exceed the number of dimensions.");
      for(size_t i=0; i<n; i++) {
        if(projectdim[i]>idim)
          mexErrMsgIdAndTxt( "MATLAB:mexController", "Cannot project controller input onto given dimensions.");
      }
      uiProj=is->project(ui,projectdim);
      size_t size=uiProj.size();
      if(size) {
        plhs[0]=mxCreateDoubleMatrix(size,(mwSize)n,mxREAL);
        double *points=mxGetPr(plhs[0]);
        for(size_t v=0;v<size;v++) {
          is->itox(uiProj[v],u);
          for(size_t w=0;w<n;w++)
            points[w*size+v]=u[projectdim[w]-1];
        }
      }
      else
        plhs[0]=mxCreateDoubleScalar(-1);
    }
    else {
      if(ui.size()) {
        plhs[0]=mxCreateDoubleMatrix(ui.size(),(mwSize)idim,mxREAL);
        double *points=mxGetPr(plhs[0]);
        for(size_t v=0;v<ui.size();v++) {
          is->itox(ui[v],u);
          for(size_t w=0;w<idim;w++)
            points[w*ui.size()+v]=u[w];
        }
      }
      else
        plhs[0]=mxCreateDoubleScalar(-1);
    }
    return;
  }
  /* number of elements in the controller domain */
  size_t dsize=game->sizeOfDomain();
  /* domain is empty */
  if(!dsize) {
    plhs[0]=mxCreateDoubleMatrix(0,0,mxREAL);
    return;
  }
  /* domain: return grid points in the domain of the controller */
  if (!strcmp(command,"domain")) {
    if(nrhs>4 && mxIsNumeric(prhs[4])) {
      mwSize n=mxGetN(prhs[4]);
      double *ptr=mxGetPr(prhs[4]);
      std::vector<size_t> projectdim;
      for(size_t i=0; i<n; i++) {
        projectdim.push_back((size_t)ptr[i]);
      }
      if(n> sdim)
        mexErrMsgIdAndTxt( "MATLAB:mexController", "Number of project dimensions cannot exceed the number of dimensions.");
      for(size_t i=0; i<n; i++) {
        if(projectdim[i]>sdim)
          mexErrMsgIdAndTxt( "MATLAB:mexController", "Cannot project domain onto given dimensions.");
      }
      std::vector<size_t> domain;
      domain=game->getDomainIdx();
      std::vector<size_t> domainProj=ss->project(domain,projectdim);
      size_t size=domainProj.size();
      /* create matrix to store points in the controller domain */
      plhs[0]=mxCreateDoubleMatrix(size,(mwSize)n,mxREAL);
      double *points=mxGetPr(plhs[0]);
      std::vector<double> x(sdim);
      for(size_t i=0; i<size; i++) {
        ss->itox(domainProj[i],x);
        for(size_t j=0; j<n; j++)
          points[i+j*size]=x[projectdim[j]-1];
      }
    } else {
      /* copy the indices of the grid points in the domain */
      size_t *domain=(size_t*)malloc(dsize*sizeof(size_t));
      game->getDomainIdx(domain);
      /* create matrix to store points in the controller domain */
      plhs[0]=mxCreateDoubleMatrix(dsize,(mwSize)sdim,mxREAL);
      double *points=mxGetPr(plhs[0]);
      std::vector<double> x(sdim);
      for(size_t i=0; i<dsize; i++) {
        ss->itox(domain[i],x);
        for(size_t j=0; j<sdim; j++)
          points[i+j*dsize]=x[j];
      }
    }
    return;
  }
  /* value: return vector containing the values of the value function at the
   * grid points in the controller domain */
  if (!strcmp(command,"value")) {
    /* create matrix to store grid points in the abstract set */
    plhs[0]=mxCreateDoubleMatrix((mwSize)(dsize),1,mxREAL);
    double *val=mxGetPr(plhs[0]);
    game->getValue(val);
    return;
  }
 

  return;
}

