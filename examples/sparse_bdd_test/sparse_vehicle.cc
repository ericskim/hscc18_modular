/*
 * vehicle.cc
 *
 *  created: Oct 2016
 *   author: Matthias Rungger
 */

/*
 * information about this example is given in
 * http://arxiv.org/abs/1503.03715
 * doi: 10.1109/TAC.2016.2593947
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
const double tau = 0.3;

/*
 * data types for the state space elements and input space
 * elements used in uniform grid and ode solvers
 */
using state_type = std::array<double,state_dim>;
using input_type = std::array<double,input_dim>;

/* abbrev of the type for abstract states and inputs */
using abs_type = scots::abs_type;

/* we integrate the vehicle ode by tau sec (the result is stored in x)  */
auto  vehicle_post = [](state_type &x, const input_type &u) {
  /* the ode describing the vehicle */
  auto rhs =[](state_type& xx,  const state_type &x, const input_type &u) {
    double alpha=std::atan(std::tan(u[1])/2.0);
    xx[0] = u[0]*std::cos(alpha+x[2])/std::cos(alpha);
    xx[1] = u[0]*std::sin(alpha+x[2])/std::cos(alpha);
    xx[2] = u[0]*std::tan(u[1]);
  };
  /* simulate (use 10 intermediate steps in the ode solver) */
  scots::runge_kutta_fixed4(rhs,x,u,state_dim,tau,10);
};

/* we integrate the growth bound by 0.3 sec (the result is stored in r)  */
auto radius_post = [](state_type &r, const state_type &, const input_type &u) {
  double c = std::abs(u[0])*std::sqrt(std::tan(u[1])*std::tan(u[1])/4.0+1);
  r[0] = r[0]+c*r[2]*tau;
  r[1] = r[1]+c*r[2]*tau;
};

int main() {
  /* to measure time */
  TicToc tt;
  /* cudd manager */
  Cudd mgr;
  mgr.AutodynEnable();
  //mgr.AutodynDisable();

  /* try to read data from files */
  scots::SymbolicSet ss_pre;
  scots::SymbolicSet ss_post;
  if(!scots::read_from_file(mgr,ss_pre,"vehicle_state_pre") ||
     !scots::read_from_file(mgr,ss_post,"vehicle_state_post")) {
   /* lower bounds of the hyper rectangle */
    state_type s_lb={{0,0,-3.5}};
    /* upper bounds of the hyper rectangle */
    state_type s_ub={{10,10,3.5}};
    /* grid node distance diameter */
    state_type s_eta={{.2,.2,.2}};
    /* construct SymbolicSet with the UniformGrid information for the state space
     * and BDD variable IDs for the pre */
    ss_pre = scots::SymbolicSet(mgr, state_dim,s_lb,s_ub,s_eta);
    /* construct SymbolicSet with the UniformGrid information for the state space
     * and BDD variable IDs for the post */
    ss_post = scots::SymbolicSet(mgr, state_dim,s_lb,s_ub,s_eta);
 
    scots::write_to_file(ss_pre,"vehicle_state_pre");
    scots::write_to_file(ss_post,"vehicle_state_post");
  }
  std::cout << "Unfiorm grid details:" << std::endl;
  ss_pre.print_info(1);

  /* try to read data from files */
  scots::SymbolicSet ss_input;
  if(!scots::read_from_file(mgr,ss_input,"vehicle_input_alphabet")) {
    /* construct grid for the input space */
    /* lower bounds of the hyper rectangle */
    input_type i_lb={{-1,-1}};
    /* upper bounds of the hyper rectangle */
    input_type i_ub={{ 1, 1}};
    /* grid node distance diameter */
    input_type i_eta={{.3,.3}};
    ss_input = scots::SymbolicSet(mgr, input_dim,i_lb,i_ub,i_eta);
    scots::write_to_file(ss_input,"vehicle_input_alphabet");
  }
  ss_input.print_info(1);

  scots::SymbolicSet set;
  /* initialize SymbolicModel class with the abstract state and input alphabet */
  scots::SymbolicModel<state_type,input_type> sym_model(ss_pre,ss_input,ss_post);
  set = scots::SymbolicSet(scots::SymbolicSet(ss_pre,ss_input),ss_post);

  /* Declare dependency and print out */
  scots::DT_Dependency dep(state_dim,input_dim);
  std::vector<std::vector<int> > state_rhs({{0,2},{1,2},{2}}), input_rhs(state_dim, {0,1});
  dep.set_rhs(state_rhs, input_rhs);
  std::cout << dep;

  /* Sparse abstraction*/
  BDD sparse_TF;
  if(!scots::read_from_file(mgr,set,sparse_TF,"vehicle_sparse_tf")) {

    std::cout << "\nComputing the sparse transition function: " << std::endl;
    tt.tic();
    size_t no_trans;
    sparse_TF = sym_model.compute_sparse_gb(mgr,vehicle_post,radius_post,dep, no_trans);
    tt.toc();
    std::cout << "Number of transitions: " << no_trans << std::endl;
    if(!getrusage(RUSAGE_SELF, &usage))
      std::cout << "Memory per transition: " << usage.ru_maxrss/(double)no_trans << std::endl;

   scots::write_to_file(mgr,set,sparse_TF,"vehicle_sparse_tf");
  }

  /*Brute force abstraction*/
  BDD TF;
  /* does there exist the transition function file ?*/
  if(!scots::read_from_file(mgr,set,TF,"vehicle_tf")) {

    std::cout << "\nComputing the transition function: " << std::endl;
    tt.tic();
    size_t no_trans;
    TF = sym_model.compute_gb(mgr,vehicle_post,radius_post,no_trans);
    tt.toc();
    std::cout << "Number of transitions: " << no_trans << std::endl;
    if(!getrusage(RUSAGE_SELF, &usage))
      std::cout << "Memory per transition: " << usage.ru_maxrss/(double)no_trans << std::endl;
    
    scots::write_to_file(mgr,set,TF,"vehicle_tf");
  }

  /* Compare */
  if(sparse_TF == TF){
    std::cout << "\nNo difference between sparse and brute force abstractions!\n"; 
  }
  else {
    std::cout << "\nSparse and brute force abstractions are different!!!\n";
  }

}
