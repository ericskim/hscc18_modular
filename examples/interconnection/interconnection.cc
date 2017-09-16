/*
 * vehicle.cc
 *
 *  created: Sep 2017
 *   author: Eric Kim
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

/* time profiling */
#include "TicToc.hh"
/* memory profiling */
#include <sys/time.h>
#include <sys/resource.h>
struct rusage usage;

/* state space dim */
const int state_dim=2;
/* input space dim */
const int control_dim=1;
/* exog space dim*/
const int exog_dim = 1;
/* input space is a cartesian product*/
const int input_dim = state_dim + control_dim + exog_dim; 
/* Create N identical systems */
const int N = 8;

/*
 * data types for the state space elements and input space
 * elements used in uniform grid and ode solvers
 */
using state_type = std::array<double,state_dim>;
using control_type = std::array<double,control_dim>;
using exog_type  = std::array<double,exog_dim>;
using input_type  = std::array<double,input_dim>;

/* abbrev of the type for abstract states and inputs */
using abs_type = scots::abs_type;

inline double saturate(double x, double lb, double ub){
  if (x < lb){
    return lb;
  } else if(x > ub){
    return ub;
  }
  else{
    return x;
  }
}



int main() {
  /* to measure time */
  TicToc tt;
  /* cudd manager */
  Cudd mgr;
  mgr.AutodynEnable();
  //mgr.AutodynDisable();

  /* Function */ 
  auto dynamics = [](state_type x, const control_type &u, const exog_type &w) -> state_type {
    state_type post;
    post[0] = saturate(x[0] + x[1], 0.0, 40.0);
    post[1] = saturate(x[1] + u[1] + w[0], -1.0, 1.0);
    return post;
  };

  /* Takes an input box and computes an overapproximating box. Both are represented with their 
  */
  auto sys_overapprox = [dynamics](const input_type &i_ll, const input_type &i_ur, state_type& o_ll, state_type& o_ur){
    o_ll = dynamics({{i_ll[0], i_ll[1]}}, {{i_ll[2]}}, {{i_ll[3]}});
    o_ur = dynamics({{i_ur[0], i_ur[1]}}, {{i_ur[2]}}, {{i_ur[3]}});
  };

  /* State spaces */
  std::vector<scots::SymbolicSet> ss_pre; ss_pre.resize(N);
  std::vector<scots::SymbolicSet> ss_post; ss_post.resize(N);
  state_type s_lb={{0, -1}};
  /* upper bounds of the hyper rectangle */
  state_type s_ub={{40, 1}};
  /* grid node distance diameter */
  state_type s_eta={{.2,.2}};
  for (int i = 0; i < N; i++){
    ss_pre[i] = scots::SymbolicSet(mgr, state_dim,s_lb,s_ub,s_eta);
    ss_post[i] = scots::SymbolicSet(mgr, state_dim,s_lb,s_ub,s_eta);
  }
  ss_pre[0].print_info(1);

  std::vector<double> x; 
  ss_pre[0].itox(0, x);
  std::cout << x[0] << " " << x[1] << "\n";
  ss_pre[0].itox(200, x);
  std::cout << x[0] << " " << x[1] << "\n";
  ss_pre[0].itox(201, x);
  std::cout << x[0] << " " << x[1] << "\n";

  /*Input spaces*/
  std::vector<scots::SymbolicSet> ss_control; ss_control.resize(N);
  /* lower bounds of the hyper rectangle */
  control_type i_lb={{-1}};
  /* upper bounds of the hyper rectangle */
  control_type i_ub={{ 1}};
  /* grid node distance diameter */
  control_type i_eta={{.2}};
  for(int i = 0; i < N; i++){
    ss_control[i] = scots::SymbolicSet(mgr, control_dim,i_lb,i_ub,i_eta);
  }


  /* Exogenous space */
  std::vector<scots::SymbolicSet> ss_exog; ss_exog.resize(N);
  exog_type e_lb = {{0}};
  exog_type e_ub = {{40}};
  exog_type e_eta = {{1.0}};
  for (int i = 0; i < N; i++){
     ss_exog[i] = scots::SymbolicSet(mgr, exog_dim,e_lb,e_ub,e_eta);
  }

  /* Declare dependencies for individual systems */
  std::vector<scots::FunctionDependency> sysdeps(N, scots::FunctionDependency());
  for (int i = 0; i < N; i++){
    sysdeps[i] = scots::FunctionDependency({ss_pre[i], ss_control[i], ss_exog[i]},{ss_post[i]});
    sysdeps[i].set_dependency(ss_post[i][0], {ss_pre[i][0],ss_pre[i][1]});
    sysdeps[i].set_dependency(ss_post[i][1], {ss_pre[i][1],ss_control[i][0],ss_exog[i][0]});
  }

  std::vector<scots::FunctionAbstracter<input_type, state_type> > abs_comp(N, scots::FunctionAbstracter<input_type, state_type>());
  std::vector<BDD> abs_systems(N, mgr.bddOne());
  /*Compute system abstractions using dependencies*/
  for (int i = 0; i < N; i++){
    abs_comp[i] = scots::FunctionAbstracter<input_type, state_type>(sysdeps[i], sys_overapprox);
    abs_systems[i] = abs_comp[i].compute_abstraction(mgr);
  }

  // std::vector<std::vector<int>> state_rhs(state_dim, std::vector<int>()), control_rhs(state_dim, std::vector<int>());
  // state_rhs[0].push_back({0});
  // state_rhs[1].push_back({1});
  // dep.set_rhs(state_rhs, control_rhs);
  // std::cout << dep;

  // scots::SymbolicSet set;

  // /* initialize SymbolicModel class with the abstract state and input alphabet */
  // scots::SymbolicModel<state_type,control_type> sym_model(ss_pre,ss_control,ss_post);

  // /** Compute transition using sparsity-aware algorithm **/
  // BDD sparse_TF;
  // if(!scots::read_from_file(mgr,set,sparse_TF,"sparse_tf")) {

  //   set = scots::SymbolicSet(scots::SymbolicSet(ss_pre,ss_control),ss_post);

  //   std::cout << "\nComputing the sparse transition function: " << std::endl;
  //   tt.tic();
  //   size_t no_trans;
  //   sparse_TF = sym_model.compute_sparse_gb(mgr,identity_post,radius_post,dep,no_trans);
  //   tt.toc();
  //   std::cout << "Number of transitions: " << no_trans << std::endl;
  //   if(!getrusage(RUSAGE_SELF, &usage))
  //     std::cout << "Memory per transition: " << usage.ru_maxrss/(double)no_trans << std::endl;
  //   //scots::write_to_file(mgr,set,sparse_TF,"sparse_tf");
  // }

  /* Compute abstraction of the averaging interconnection */

}

