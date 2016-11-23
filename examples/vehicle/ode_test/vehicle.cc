/****************************************************************************
	Driver file of the dp_tides program
	This file has been created by MathTIDES (2.00) October 28, 2015, 11:44

	Copyright (C) 2010 A. Abad, R. Barrio, F. Blesa, M. Rodriguez
	Grupo de Mecanica Espacial
	University of Zaragoza
	SPAIN

	http://gme.unizar.es/software/tides
	Contact: <tides@unizar.es>

	This file is part of TIDES.

	TIDES is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	TIDES is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with TIDES.  If not, see <http://www.gnu.org/licenses/>.

*****************************************************************************/
#include <iostream>
#include <array>

#include "UniformGrid.hh"
#include "Abstraction.hh"

#include "TicToc.hh"


extern "C" {
  #include "dp_tides.h"
  #include "vehicleODE.h"
}


/* state space dim */
#define sDIM 3
#define iDIM 2

/* functions to setup the state space and input space of the vehicle example */

/* data types for the ode solver */
typedef std::array<double,3> state_type;
typedef std::array<double,2> input_type;

/* we integrate the vehicle ode by 0.3 sec (the result is stored in x)  */
auto  vehicle_post_hi = [](state_type &x, input_type &u) -> void {
	double tolrel = 1.e-16 ;
	double tolabs = 1.e-16;

  int nipt = 1;

	double dt = 0.3; /* (tend-tini)/nipt */

	dp_tides_delta(vehicleODE, NULL, sDIM, iDIM, 0, &x[0], &u[0], 0, dt, nipt, tolrel, tolabs, NULL, NULL);

};

/* we integrate the vehicle ode by 0.3 sec (the result is stored in x)  */
auto  vehicle_post = [](state_type &x, input_type &u) -> void {

  /* the ode describing the vehicle */
  auto rhs =[](state_type& xx,  const state_type &x, input_type &u) -> void {
      double alpha=std::atan(std::tan(u[1])/2.0);
      xx[0] = u[0]*std::cos(alpha+x[2])/std::cos(alpha);
      xx[1] = u[0]*std::sin(alpha+x[2])/std::cos(alpha);
      xx[2] = u[0]*std::tan(u[1]);
  };


  /* runge kutte order 4 */
  state_type k[4];
  state_type tmp;

  size_t nint=3; /* number of intermediate step size */
  double h=0.1; /* h* nint = sampling time */

  for(size_t t=0; t<nint; t++) {
    rhs(k[0],x, u);
    for(size_t i=0;i<sDIM;i++)
      tmp[i]=x[i]+h/2*k[0][i];

    rhs(k[1],tmp, u);
    for(size_t i=0;i<sDIM;i++)
      tmp[i]=x[i]+h/2*k[1][i];

    rhs(k[2],tmp, u);
    for(size_t i=0;i<sDIM;i++)
      tmp[i]=x[i]+h*k[2][i];

    rhs(k[3],tmp, u);
    for(size_t i=0; i<sDIM; i++)
      x[i] = x[i] + (h/6)*(k[0][i] + 2*k[1][i] + 2*k[2][i] + k[3][i]);
  }
};

int main() {
  /* to measure time */
  TicToc tt;

  /****************************************************************************/
  /* construct grid for the state space */
  /****************************************************************************/
  /* setup the workspace of the synthesis problem and the uniform grid */
  /* lower bounds of the hyper rectangle */
  double lb[sDIM]={0,0,-M_PI-0.4};  
  /* upper bounds of the hyper rectangle */
  double ub[sDIM]={10,10,M_PI+0.4}; 
  /* grid node distance diameter */
  double eta[sDIM]={.2,.2,.2};   
  scots::UniformGrid ss(sDIM,lb,ub,eta);
  std::cout << "Unfiorm grid details:" << std::endl;

  /****************************************************************************/
  /* construct grid for the input space */
  /****************************************************************************/
  /* lower bounds of the hyper rectangle */
  double ilb[sDIM]={-1,-1};  
  /* upper bounds of the hyper rectangle */
  double iub[sDIM]={1,1}; 
  /* grid node distance diameter */
  double ieta[sDIM]={.3,.3};   
  scots::UniformGrid is(iDIM,ilb,iub,ieta);

 

  size_t N=ss.getN();
  size_t M=is.getN();

   /* state and input variables */
  state_type x;
  state_type y;
  input_type u;

 
  tt.tic();

  /* loop over all states */
  for(size_t i=0; i<N; i+=10) {
    for(size_t j=0; j<M; j++) {
      /* current state */
      ss.idxToElement(i,&x[0]);
      ss.idxToElement(i,&y[0]);
      /* current input */
      is.idxToElement(j,&u[0]);

      vehicle_post_hi(x,u);
      vehicle_post(y,u);

      for(size_t k=0; k<sDIM; k++) {
        
        if(std::abs(y[k]-x[k])>1e-8) {
          std::cout << x[0] << " " << x[1] << " " << x[2] << std::endl;
          std::cout << y[0] << " " << y[1] << " " << y[2] << std::endl;
        }

      }

     
 
    }
  }


  tt.toc(); 

  return 1;
}

