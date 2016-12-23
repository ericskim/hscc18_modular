/*
 * aircraft.cc
 *
 *  created on: 18.01.2016
 *  author: rungger
 *
 */

/*
 * information about this example is given in
 * http://arxiv.org/abs/1503.03715
 * doi: 10.1109/TAC.2016.2593947
 * doi: 10.1109/CDC.2015.7403185
 *
 */

#include <iostream>
#include <array>

/* SCOTS header */
#include "UniformGrid.hh"
#include "TransitionSystem.hh"
#include "AbstractionGB.hh"
#include "ReachabilityGame.hh"

/* time profiling */
#include "TicToc.hh"

/* ode solver */
#include "RungeKutta4.hh"

/* ode solver */
OdeSolver ode_solver;

/* state space dim */
const int sDIM=3;
const int iDIM=2;

/* sampling time */
const double tau = 0.25;

/*
 * data types of the state space elements and input 
 * space elements used in uniform grid and ode solver
 *
 */
typedef std::array<double,sDIM> state_type;
typedef std::array<double,iDIM> input_type;

/* we integrate the aircraft ode by 0.25 sec (the result is stored in x)  */
double mg = 60000.0*9.81;
double mi = 1.0/60000;
auto aircraft_post = [] (state_type &x, const input_type &u) {
  /* the ode describing the aircraft */
  auto rhs =[] (state_type& xx,  const state_type &x, const input_type &u) {
    double c=(1.25+4.2*u[1]);
    xx[0] = mi*(u[0]*std::cos(u[1])-(2.7+3.08*c*c)*x[0]*x[0]-mg*std::sin(x[1]));
    xx[1] = (1.0/(60000*x[0]))*(u[0]*std::sin(u[1])+68.6*c*x[0]*x[0]-mg*std::cos(x[1]));
    xx[2] = x[0]*std::sin(x[1]);
  };
  /* use 10 intermediate steps (check ./helper/ode_test to find parameters) */
  ode_solver(rhs,x,u,sDIM,tau,10);
};

double L[3][2];
/* to account for input disturbances */
state_type w={{.108,0.002,0}};
/* we integrate the growth bound by 0.25 sec (the result is stored in r)  */
auto radius_post = [] (state_type &r, const state_type &, const input_type &u) {

  /* lipschitz matrix computed with mupad/mathematica check the ./helper directory */
  L[0][0]=-0.001919*(2.7+3.08*(1.25+4.2*u[1])*(1.25+4.2*u[1]));
  L[0][1]=9.8100000000000041;

  L[1][0]=0.00294+0.00481*u[1];
  L[1][1]=0.003526;

  L[2][0]=0.074252;
  L[2][1]=83.2;

  /* the ode for the growth bound */
  auto rhs =[] (state_type& rr,  const state_type &r, const input_type &) {
    rr[0] = L[0][0]*r[0]+L[0][1]*r[1]+w[0]; /* L[0][2]=0 */
    rr[1] = L[1][0]*r[0]+L[1][1]*r[1]+w[1]; /* L[1][2]=0 */
    rr[2] = L[2][0]*r[0]+L[2][1]*r[1]+w[2]; /* L[2][2]=0 */
  };
  /* use 10 intermediate steps (check ./helper/ode_test to find parameters) */
  ode_solver(rhs,r,u,sDIM,tau,10);
};

int main() {
  /* to measure time */
  TicToc tt;

  /****************************************************************************/
  /* construct grid for the state space */
  /****************************************************************************/
  /* setup the workspace of the synthesis problem and the uniform grid */
  /* grid node distance diameter */
  /* optimized values computed according to doi: 10.1109/CDC.2015.7403185 */
  state_type eta={{25.0/362,3*M_PI/180/66,56.0/334}}; 
  /* lower bounds of the hyper rectangle */
  state_type lb={{58,-3*M_PI/180,0}};
  /* upper bounds of the hyper rectangle */
  state_type ub={{83,0,56}}; 
  /* measurement disturbances  */
  state_type z={{0.0125,0.0025/180*M_PI,0.05}};
  scots::UniformGrid<state_type> ss(sDIM,lb,ub,eta,z);
  std::cout << "Unfiorm grid details:" << std::endl;
  ss.printInfo(1);

  /****************************************************************************/
  /* construct grid for the input space */
  /****************************************************************************/
  /* lower bounds of the hyper rectangle */
  input_type ilb={{0,0}};
  /* upper bounds of the hyper rectangle */
  input_type iub={{32000,8*M_PI/180}};
  /* grid node distance diameter */
  input_type ieta={{32000,8.0/9.0*M_PI/180}};

  scots::UniformGrid<input_type> is(iDIM,ilb,iub,ieta);
  is.printInfo(1);

  /* transition system to be computed */
  scots::TransitionSystem ts;


  tt.tic();
  scots::AbstractionGB<state_type,input_type> abs(ss,is,ts);

  abs.computeTransitionRelation(aircraft_post, radius_post);

  std::cout << "Number of transitions: " << ts.getNoTransitions() << std::endl;
  tt.toc();
  state_type x;
  auto target = [&](const size_t idx)->bool {
    ss.itox(idx,x);
    /* function returns 1 if cell associated with x is in target set  */
    if(  63 <= (x[0]-eta[0]/2.0) && (x[0]+eta[0]/2.0)<= 75 &&
         -3*M_PI/180 <= (x[1]-eta[1]/2.0) &&  (x[1]+eta[1]/2.0) <= 0 &&
         0 <= (x[2]-eta[2]/2.0) &&  (x[2]+eta[2]/2.0) <= 2.5 &&
         -0.91 <=  ( (x[0]+eta[0]/2.0) * std::sin(x[1]-eta[1]/2.0) )
      )
      return true;
    return false;
  };

  tt.tic();

  std::cout << "Solve game " << std::endl;

  scots::ReachabilityGame reach(ts);
  reach.solve(target);
  tt.toc();

  std::cout << "Size: " << reach.size() << std::endl;


  return 1;
}
