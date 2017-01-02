/*
 * vehicle.cc
 *
 *  created on: 26.10.2015
 *      author: rungger
 */

/*
 * information about this example is given in
 * http://arxiv.org/abs/1503.03715
 * doi: 10.1109/TAC.2016.2593947
 *
 */

#include <iostream>
#include <array>
#include <iomanip>

/* SCOTS header */
#include "UniformGrid.hh"
#include "TransitionSystem.hh"
#include "AbstractionGB.hh"
#include "ReachabilityGame.hh"
#include "StaticController.hh"

/* time profiling */
#include "TicToc.hh"

/* ode solver */
#include "RungeKutta4.hh"
OdeSolver ode_solver;

/* state space dim */
const int state_dim=3;
/* input space dim */
const int input_dim=2;

/* sampling time */
const double tau = 0.3;


/*
 * data types for the state space elements and input space
 * elements used in uniform grid and ode solvers
 *
 */
typedef std::array<double,state_dim> state_type;
typedef std::array<double,input_dim> input_type;

/* we integrate the vehicle ode by 0.3 sec (the result is stored in x)  */
auto  vehicle_post = [](state_type &x, const input_type &u) {
  /* the ode describing the vehicle */
  auto rhs =[](state_type& xx,  const state_type &x, const input_type &u) {
    double alpha=std::atan(std::tan(u[1])/2.0);
    xx[0] = u[0]*std::cos(alpha+x[2])/std::cos(alpha);
    xx[1] = u[0]*std::sin(alpha+x[2])/std::cos(alpha);
    xx[2] = u[0]*std::tan(u[1]);
  };
  /* simulate (use 10 intermediate steps in the ode solver) */
  ode_solver(rhs,x,u,state_dim,tau,10);
};

/* we integrate the growth bound by 0.3 sec (the result is stored in r)  */
auto radius_post = [](state_type &r, const state_type &, const input_type &u) {
  double c = std::abs(u[0])*std::sqrt(std::tan(u[1])*std::tan(u[1])/4.0+1);
  r[0] = r[0]+c*r[2]*tau+1e-16;
  r[1] = r[1]+c*r[2]*tau+1e-16;
};

int main() {

  /* to measure time */
  TicToc tt;

  /****************************************************************************/
  /* construct grid for the state space */
  /****************************************************************************/
  /* setup the workspace of the synthesis problem and the uniform grid */
  /* lower bounds of the hyper rectangle */
  state_type lb={{0,0,-M_PI-0.4}};
  /* upper bounds of the hyper rectangle */
  state_type ub={{10,10,M_PI+0.4}};
  /* grid node distance diameter */
  state_type eta={{.1,.1,.2}};
  scots::UniformGrid<state_type> ss(state_dim,lb,ub,eta);
  std::cout << "Unfiorm grid details:" << std::endl;
  ss.printInfo(1);

  /****************************************************************************/
  /* construct grid for the input space */
  /****************************************************************************/
  /* lower bounds of the hyper rectangle */
  input_type ilb={{-1,-1}};
  /* upper bounds of the hyper rectangle */
  input_type iub={{1,1}};
  /* grid node distance diameter */
  input_type ieta={{.3,.3}};
  scots::UniformGrid<input_type> is(input_dim,ilb,iub,ieta);

  /****************************************************************************/
  /* set up constraint functions with obtacles */
  /****************************************************************************/
  double H[15][4] = {
    { 1  , 1.2, 0  ,   9 },
    { 2.2, 2.4, 0  ,   5 },
    { 2.2, 2.4, 6  ,  10 },
    { 3.4, 3.6, 0  ,   9 },
    { 4.6, 4.8, 1  ,  10 },
    { 5.8, 6  , 0  ,   6 },
    { 5.8, 6  , 7  ,  10 },
    { 7  , 7.2, 1  ,  10 },
    { 8.2, 8.4, 0  ,  8.5},
    { 8.4, 9.3, 8.3,  8.5},
    { 9.3, 10 , 7.1,  7.3},
    { 8.4, 9.3, 5.9,  6.1},
    { 9.3, 10 , 4.7,  4.9},
    { 8.4, 9.3, 3.5,  3.7},
    { 9.3, 10 , 2.3,  2.5}
  };

  /* overflow function returns 1 if x \in overflow symbol  */
  state_type x;
  //auto overflow = [&](const size_t idx) {
  auto overflow = [&](const state_type &x) {
    //ss.itox(idx,x);
    double c1= eta[0]/2.0;
    double c2= eta[1]/2.0;
    for(size_t i=0; i<15; i++) {
      if ((H[i][0]-c1) <= x[0] && x[0]<= (H[i][1]+c1) && (H[i][2]-c2) <= x[1] && x[1] <= (H[i][3]+c2))
        return true;
    }
    return false;
  };

  /* transition system to be computed */
  scots::TransitionSystem ts;

  tt.tic();
  scots::AbstractionGB<state_type,input_type> abs(ss,is,ts);

  abs.computeTransitionRelation(vehicle_post, radius_post, overflow);

  std::cout << "Number of transitions: " << ts.getNoTransitions() << std::endl;
  tt.toc();

  /* define target set */
  auto target = [&](size_t idx) {
    ss.itox(idx,x);
    /* function returns 1 if cell associated with x is in target set  */
    if (9 <= (x[0]-eta[0]/2.0) && (x[0]+eta[0]/2.0)<= 9.5 && 0 <= (x[1]-eta[1]/2.0) &&  (x[1]+eta[1]/2.0) <= 0.5)
      return true;
    else
    return false;
  };

  scots::ReachabilityGame reach(ts);
  tt.tic();
  std::cout << "Solve game " << std::endl;
  reach.solve(target);
  tt.toc();
  std::cout << "Size: " << reach.size() << std::endl;

  return 1;
}
