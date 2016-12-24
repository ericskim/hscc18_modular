/*
 * test.cc
 *
 *  created on: 20.01.2016
 *      author: rungger
 */

/*
 * a computation to compare the RK4 ode solver with validated solution solver vnodelp
 *
 */
#include <iostream>
#include <array>
#include <cmath>
/* SCOTS header */
#include "UniformGrid.hh"
#include "TransitionSystem.hh"
#include "AbstractionGB.hh"
#include "ReachabilityGame.hh"
/* time profiling */
#include "TicToc.hh"
/* ode solver */
#include "RungeKutta4.hh"
/* vnodelp solver */
#include "vnode.h"

/* ode solver */
OdeSolver ode_solver;

/* state space dim */
const int state_dim=3;
/* input space dim */
const int input_dim=2;
/* sampling time */
const double tau = 0.25;

/* setup Runge Kutta solver */
typedef std::array<double,state_dim> state_type;
typedef std::array<double,input_dim> input_type;

/* setup rhs for  vnodelp solver */
auto rhs =[] (state_type& xx,  const state_type &x, const input_type &u) {

  double mg = 60000.0*9.81;
  double mi = 1.0/60000;
  double c=(1.25+4.2*u[1]);

  xx[0] = mi*(u[0]*std::cos(u[1])-(2.7+3.08*c*c)*x[0]*x[0]-mg*std::sin(x[1]));
  xx[1] = (1.0/(60000*x[0]))*(u[0]*std::sin(u[1])+68.6*c*x[0]*x[0]-mg*std::cos(x[1]));
  xx[2] = x[0]*std::sin(x[1]);
};



/* setup rhs for  vnodelp solver */
template<typename var_type> 
void aircraft(int n, var_type* xx, const var_type* x, var_type t, void* param) {

  /* inputs are passed by parameters */
  vnodelp::iVector u(input_dim);
  u[0] = ((vnodelp::interval*)param)[0];
  u[1] = ((vnodelp::interval*)param)[1];

  vnodelp::interval mg = vnodelp::interval(60000.0)*9.81;
  vnodelp::interval mi = vnodelp::interval(1.0)/60000.0;
  vnodelp::interval c = vnodelp::interval(1.25+4.2*u[1]);

  xx[0] = mi*(u[0]*cos(u[1])-(2.7+3.08*c*c)*x[0]*x[0]-mg*sin(x[1]));
  xx[1] = (1.0/(60000.0*x[0]))*(u[0]*sin(u[1])+68.6*c*x[0]*x[0]-mg*cos(x[1]));
  xx[2] = x[0]*sin(x[1]);

};

int main() {
  /* to measure time */
  TicToc tt;

  /* define accuracy you want to check the Runge Kutta solver for */
  double accuracy =  1e-10;

  /****************************************************************************/
  /* construct grid points to check for ode solution */
  /****************************************************************************/
  /* grid node distance diameter */
  state_type eta={{5.0/42,M_PI/12600,4.0/15}};   
  /* lower bounds of the hyper rectangle */
  state_type lb={{58,-3*M_PI/180,0}};  
  /* upper bounds of the hyper rectangle */
  state_type ub={{83,0,56}}; 
  scots::UniformGrid<state_type> ss(state_dim,lb,ub,eta);
  std::cout << "Unfiorm grid details:" << std::endl;
  ss.printInfo(1);

  /****************************************************************************/
  /* construct grid for the input space */
  /****************************************************************************/
  /* lower bounds of the hyper rectangle */
  input_type ilb={{0,0}};  
  /* upper bounds of the hyper rectangle */
  input_type iub={{36000,8*M_PI/180}};  
  /* grid node distance diameter */
  input_type ieta={{36000,8.0/9.0*M_PI/180}};  
  scots::UniformGrid<input_type> is(input_dim,ilb,iub,ieta);

  /* variables used in the integration */
  state_type x;
  input_type u;

  vnodelp::iVector ix(state_dim);
  vnodelp::iVector iu(state_dim);

  /* set up vnodelp solver */
  vnodelp::AD* ad = new vnodelp::FADBAD_AD(state_dim,aircraft,aircraft,(void*)&iu);
  vnodelp::VNODE* vnode_solver = new vnodelp::VNODE(ad);
  vnode_solver->setTols(1e-16,1e-16);

  size_t N=ss.getN();
  size_t M=is.getN();

  tt.tic();

  /* loop over all states and inputs */
  for(size_t i=0; i<N/1000; i++) {
    for(size_t j=0; j<M; j++) {
      
      /* reset internal data of ode solver */
      vnode_solver->setFirstEntry();
      vnodelp::interval t0 = 0.0;
      vnodelp::interval itau = tau;

      /* get current grid point */
      ss.itox(i,x);
      for(size_t k=0; k<state_dim; k++) 
        ix[k] = x[k];

      /* current input */
      is.itox(j,u);
      iu[0]=u[0];
      iu[1]=u[1];

      /* tell solver that params have changed */
      ad->eval((void*)&iu); 

      /* solve ode with vnodelp solver */
      vnode_solver->integrate(t0,ix,itau);
      if(!vnode_solver->successful())
        std::cout<<"VNODE-LP something wrong " << std::endl;

      /* solve ode with RungeKutta solver */
      ode_solver(rhs,x,u,state_dim,tau,10);

      /* check if RungeKutta accuracy exceeds the defined accuracy */
      for(size_t k=0; k<state_dim; k++) {
        double diff_left = std::abs(x[k]-vnodelp::inf(ix[k]));
        double diff_right = std::abs(x[k]-vnodelp::sup(ix[k]));
        
        if((diff_left >= accuracy || diff_right>= accuracy) && vnodelp::width(ix[k])>accuracy) {
          std::cout<<"RungeKutta solver error is not guaranteed  to be within " << 2*accuracy << std::endl;
          delete vnode_solver;
          delete ad;
          return 0;
        }
      }
    }
  }

  std::cout<<"RungeKutta solver error is within " << 2*accuracy << std::endl;

  tt.toc(); 

  delete vnode_solver;
  delete ad;

  return 1;
}


