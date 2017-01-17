/*
 * aircraft.cc
 *
 *  created: Oct 2016
 *  author: Matthias Rungger
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
const int state_dim=3;
/* input space dim */
const int input_dim=2;
/* sampling time */
const double tau = 0.25;

/* data types of the state space elements and input 
 * space elements used in uniform grid and ode solver */
using state_type = std::array<double,state_dim>;
using input_type = std::array<double,input_dim>;

/* we integrate the aircraft ode by 0.25 sec (the result is stored in x)  */
auto aircraft_post = [] (state_type &x, const input_type &u) {
  /* the ode describing the aircraft */
  auto rhs =[] (state_type& xx,  const state_type &x, const input_type &u) {
    double mg = 60000.0*9.81;
    double mi = 1.0/60000;
    double c=(1.25+4.2*u[1]);
    xx[0] = mi*(u[0]*std::cos(u[1])-(2.7+3.08*c*c)*x[0]*x[0]-mg*std::sin(x[1]));
    xx[1] = (1.0/(60000*x[0]))*(u[0]*std::sin(u[1])+68.6*c*x[0]*x[0]-mg*std::cos(x[1]));
    xx[2] = x[0]*std::sin(x[1]);
  };
  /* use 10 intermediate steps */
  scots::runge_kutta_fixed4(rhs,x,u,state_dim,tau,10);
};

/* we integrate the growth bound by 0.25 sec (the result is stored in r)  */
auto radius_post = [] (state_type &r, const state_type &, const input_type &u) {
  /* lipschitz matrix computed with mupad/mathematica check the ./helper directory */
  double L[3][2];
  L[0][0]=-0.00191867*(2.7+3.08*(1.25+4.2*u[1])*(1.25+4.2*u[1]));
  L[0][1]=9.81;
  L[1][0]=0.002933+0.004802*u[1];
  L[1][1]=0.003623;
  L[2][0]=0.07483;
  L[2][1]=83.22;
  /* to account for input disturbances */
  state_type w={{.108,0.002,0}};
  /* the ode for the growth bound */
  auto rhs =[&] (state_type& rr,  const state_type &r, const input_type &) {
    rr[0] = L[0][0]*r[0]+L[0][1]*r[1]+w[0]; /* L[0][2]=0 */
    rr[1] = L[1][0]*r[0]+L[1][1]*r[1]+w[1]; /* L[1][2]=0 */
    rr[2] = L[2][0]*r[0]+L[2][1]*r[1]+w[2]; /* L[2][2]=0 */
  };
  /* use 10 intermediate steps */
  scots::runge_kutta_fixed4(rhs,r,u,state_dim,tau,10);
};

int main() {
  /* to measure time */
  TicToc tt;

  /* construct grid for the state space */
  /* setup the workspace of the synthesis problem and the uniform grid */
  /* grid node distance diameter */
  /* optimized values computed according to doi: 10.1109/CDC.2015.7403185 */
  state_type s_eta={{25.0/362,3*M_PI/180/66,56.0/334}}; 
  /* lower bounds of the hyper rectangle */
  state_type s_lb={{58,-3*M_PI/180,0}};
  /* upper bounds of the hyper rectangle */
  state_type s_ub={{83,0,56}}; 
  scots::UniformGrid ss(state_dim,s_lb,s_ub,s_eta);
  std::cout << "Unfiorm grid details:" << std::endl;
  ss.print_info();

  /* construct grid for the input space */
  /* lower bounds of the hyper rectangle */
  input_type i_lb={{0,0}};
  /* upper bounds of the hyper rectangle */
  input_type i_ub={{32000,8*M_PI/180}};
  /* grid node distance diameter */
  input_type i_eta={{32000,8.0/9.0*M_PI/180}};
  scots::UniformGrid is(input_dim,i_lb,i_ub,i_eta);
  is.print_info();

  /* transition function of symbolic model */
  scots::TransitionFunction tf;

  /* setup object to compute the transition function */
  scots::AbstractionGB<state_type,input_type> abs(ss,is);
  /* measurement disturbances  */
  state_type z={{0.0125,0.0025/180*M_PI,0.05}};
  abs.set_measurement_error_bound(z);

  std::cout << "Computing the transition function: " << std::endl;
  tt.tic();
  abs.compute(tf,aircraft_post,radius_post);
  tt.toc();
  if(!getrusage(RUSAGE_SELF, &usage))
    std::cout << "Memory per transition: " << usage.ru_maxrss/(double)tf.get_no_transitions() << std::endl;
  std::cout << "Number of transitions: " << tf.get_no_transitions() << std::endl;

  tt.tic();
  write_to_file(tf,"tf.scs");
  tt.toc();

  /* define target set */
  state_type x;
  auto target = [&](const scots::abs_type idx) {
    ss.itox(idx,x);
    /* function returns 1 if cell associated with x is in target set  */
    if(         63 <= (x[0]-s_eta[0]/2.0) &&  (x[0]+s_eta[0]/2.0) <=  75 &&
       -3*M_PI/180 <= (x[1]-s_eta[1]/2.0) &&  (x[1]+s_eta[1]/2.0) <=   0 &&
                 0 <= (x[2]-s_eta[2]/2.0) &&  (x[2]+s_eta[2]/2.0) <= 2.5 &&
             -0.91 <= ((x[0]+s_eta[0]/2.0) * std::sin(x[1]-s_eta[1]/2.0) )) {
      return true;
    }
    return false;
  };

 
  std::cout << "\nSynthesis: " << std::endl;
  tt.tic();
  scots::WinningDomain win=scots::solve_reachability_game(tf,target);
  tt.toc();
  std::cout << "Winning domain size: " << win.get_size() << std::endl;

  std::cout << "\nWrite controller to controller.scs \n";
  if(write_to_file(scots::StaticController(ss,is,std::move(win)),"controller.scs"))
    std::cout << "Done. \n";

  return 1;
}
