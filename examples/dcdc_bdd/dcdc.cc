/*
 * dcdc.cc
 *
 *  created: Jan 2017
 *   author: Matthias Rungger
 */

/*
 * information about this example is given in the readme file
 */

#include <iostream>
#include <array>
#include <cmath>

/* SCOTS header */
#include "scots.hh"
/* ode solver */
#include "RungeKutta4.hh"


/* time profiling */
#include "TicToc.hh"
/* memory profiling */
#include <sys/time.h>
#include <sys/resource.h>
struct rusage usage;


/* state space dim */
const int state_dim=2;
/* input space dim */
const int input_dim=1;
/* sampling time */
const double tau = 0.5;

/*
 * data types for the elements of the state space 
 * and input space used by the ODE solver
 */
using state_type = std::array<double,state_dim>;
using input_type = std::array<double,input_dim>;

/* abbrev of the type for abstract states and inputs */
using abs_type = scots::abs_type;

/* parameters for system dynamics */
const double xc=70;
const double xl=3;
const double rc=0.005;
const double rl=0.05;
const double ro=1;
const double vs=1;
/* parameters for radius calculation */
const double mu=std::sqrt(2);
/* we integrate the dcdc ode by 0.5 sec (the result is stored in x)  */
auto system_post = [](state_type &x, const input_type &u) noexcept {
  /* the ode describing the dcdc converter */
  auto rhs =[](state_type& xx,  const state_type &x, const input_type &u) noexcept {
    if(u[0]==1) {
      xx[0]=-rl/xl*x[0]+vs/xl;
      xx[1]=-1/(xc*(ro+rc))*x[1];
    } else {
      xx[0]=-(1/xl)*(rl+ro*rc/(ro+rc))*x[0]-(1/xl)*ro/(5*(ro+rc))*x[1]+vs/xl;
      xx[1]=(1/xc)*5*ro/(ro+rc)*x[0]-(1/xc)*(1/(ro+rc))*x[1];
    }
	};
  scots::runge_kutta_fixed4(rhs,x,u,state_dim,tau,5);
};
/* we integrate the growth bound by 0.5 sec (the result is stored in r)  */
auto radius_post = [](state_type &r, const state_type&, const input_type &u) noexcept {
  /* the ode for the growth bound */
  auto rhs =[](state_type& rr,  const state_type &r, const input_type &u) noexcept {
    if(u[0]==1) {
      rr[0]=-rl/xl*r[0];
      rr[1]=-1/(xc*(ro+rc))*r[1];
    } else {
      rr[0]=-(1/xl)*(rl+ro*rc/(ro+rc))*r[0]+(1/xl)*ro/(5*(ro+rc))*r[1];
      rr[1]=5*(1/xc)*ro/(ro+rc)*r[0]-(1/xc)*(1/(ro+rc))*r[1];
    }
	};
  scots::runge_kutta_fixed4(rhs,r,u,state_dim,tau,5);
};

int main() {
  /* to measure time */
  TicToc tt;
  /* BDD manager */
  Cudd manager;

  manager.AutodynEnable();
  
  /* setup the workspace of the synthesis problem and the uniform grid */
  /* grid node distance diameter */
  state_type eta={{2.0/4e3,2.0/4e3}};
  /* lower bounds of the hyper-rectangle */
  state_type lb={{1.15-eta[0]/2,5.45-eta[1]/2}};
  /* upper bounds of the hyper-rectangle */
  state_type ub={{1.55+eta[0]/2,5.85+eta[1]/2}};
  scots::SymbolicSet ss_pre(manager,state_dim,lb,ub,eta);
  scots::SymbolicSet ss_post(manager,state_dim,lb,ub,eta);
  std::cout << "Unfiorm grid details:\n";
  ss_pre.print_info();

  /* construct grid for the input alphabet */
  /* hyper-rectangle [1,2] with grid node distance 1 */
  scots::SymbolicSet ss_input(manager,input_dim,input_type{{.99}},input_type{{2.1}},input_type{{1}});
  ss_input.print_info();
  
  /* compute transition function of symbolic model */
  std::cout << "Computing the transition function:\n";
  /* SymbolicModel class to compute the BDD encoding the transition function */ 
  scots::SymbolicModel<state_type,input_type> sym_model(ss_pre,ss_input,ss_post);

  tt.tic();
  size_t no_trans;
  BDD tf = sym_model.compute_gb(manager,system_post,radius_post,no_trans);
  tt.toc();

  std::cout << "No of Transitions " << no_trans  << "\n";
  if(!getrusage(RUSAGE_SELF, &usage)) {
    std::cout << "Memory pro Transition: " << usage.ru_maxrss/(double)no_trans<< "\n";
  }

  return 1;
}

