/*
 * runningmax.cc
 *
 *  created: Sep 2017
 *   author: Eric Kim
 */

/*
* Computes a max of N numbers through a running maximum, similar to the for loop below by in a declarative setting
* z = x[0]; 
* for(i = 1; i< N; i++) 
*   z = max(z, x[i]);
*
*
 */

#include <iostream>
#include <array>
#include <cmath>
#include <algorithm>

/* SCOTS header */
#include "scots.hh"

/* time profiling */
#include "TicToc.hh"
/* memory profiling */
#include <sys/time.h>
#include <sys/resource.h>
struct rusage usage;

/* state space dim */
const int state_dim=1;
/* input space dim */
const int control_dim=1;
/* exog space dim*/
const int exog_dim = 1;
/* input space of system is a cartesian product*/
const int input_dim = state_dim + control_dim + exog_dim; 
/* Create N identical systems */
const int N = 10;

const int inter_dim = N - 2; 

/*
 * data types for the state space elements and input space
 * elements used ll uniform grid
 */
using state_type = std::array<double,state_dim>;
using control_type = std::array<double,control_dim>;
using exog_type  = std::array<double,exog_dim>;
using input_type  = std::array<double,input_dim>;
using inter_type = std::array<double, inter_dim>;

/* Data types for the interconnection relation */
using prod_state_type = std::array<double, state_dim * N>;
using prod_exog_type = std::array<double, exog_dim * N>;

/* abbrev of the type for abstract states and inputs */
using abs_type = scots::abs_type;

/* Saturation function */
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


void print_support(const Cudd& mgr, const BDD& x){
  std::vector< unsigned int >  indices = mgr.SupportIndices({x});
  for (size_t i = 0; i < indices.size(); i++){
    std::cout << indices[i] << " ";
  }
  std::cout << std::endl;
}

/* Generate max functions */
void max2(std::array<double,2> ll, std::array<double, 2> ur, std::array<double, 1> &o_ll, std::array<double,1> &o_ur){
  o_ll[0] = std::max(ll[0], ll[1]); 
  o_ur[0] = std::max(ur[0], ur[1]);
}

int main() {
  /* to measure time */
  TicToc tt;
  /* cudd manager */
  Cudd mgr;
  mgr.AutodynEnable(CUDD_REORDER_SIFT_CONVERGE);
  // mgr.AutodynEnable(CUDD_REORDER_RANDOM_PIVOT);
  mgr.EnableReorderingReporting();
  //mgr.AutodynDisable();

  /* Dynamics for individual subsystem */ 
  auto dynamics = [](const state_type x, const control_type u, const exog_type w) -> state_type {
    state_type post;
    post[0] = std::min(.75*(x[0] + u[0]), w[0] + 1);
    post[0] = saturate(post[0], 0, 31.0);
    return post;
  };

  /* Takes an input box and computes an overapproximating box. Both are represented with their lower left
  and upper right corners.
  */
  auto sys_overapprox = [dynamics](const input_type i_ll, const input_type i_ur, state_type& o_ll, state_type& o_ur){
    o_ll = dynamics({i_ll[0]}, {i_ll[1]}, {i_ll[2]});
    o_ur = dynamics({i_ur[0]}, {i_ur[1]}, {i_ur[2]});
  };

  /* State spaces */
  std::vector<scots::SymbolicSet> ss_pre; ss_pre.resize(N);
  scots::SymbolicSet pre_product = scots::SymbolicSet();
  std::vector<scots::SymbolicSet> ss_post; ss_post.resize(N);
  scots::SymbolicSet post_product = scots::SymbolicSet();
  state_type s_lb={{0}};
  /* upper bounds of the hyper rectangle */
  state_type s_ub={{31}};
  /* grid node distance diameter */
  state_type s_eta={{1.0}};
  for (int i = 0; i < N; i++){
    ss_pre[i] = scots::SymbolicSet(mgr, state_dim,s_lb,s_ub,s_eta);
    ss_post[i] = scots::SymbolicSet(mgr, state_dim,s_lb,s_ub,s_eta);

    pre_product = scots::SymbolicSet(pre_product, ss_pre[i]);
    post_product = scots::SymbolicSet(post_product, ss_post[i]);
  }
  std::cout << "Pre State Product Information" << std::endl;
  pre_product.print_info(1);
  
  /*Input spaces*/
  std::vector<scots::SymbolicSet> ss_control; ss_control.resize(N);
  scots::SymbolicSet control_product = scots::SymbolicSet();
  /* lower bounds of the hyper rectangle */
  control_type i_lb={{0}};
  /* upper bounds of the hyper rectangle */
  control_type i_ub={{7}};
  /* grid node distance diameter */
  control_type i_eta={{1.0}};
  for(int i = 0; i < N; i++){
    ss_control[i] = scots::SymbolicSet(mgr, control_dim,i_lb,i_ub,i_eta);
    control_product = scots::SymbolicSet(control_product, ss_control[i]);
  }
  std::cout << "Control Product Information" << std::endl;
  control_product.print_info(1);

  /* Exogenous spaces */
  scots::SymbolicSet ss_exog;
  exog_type e_lb = {{0}};
  exog_type e_ub = {{31}};
  exog_type e_eta = {{1}};
  ss_exog = scots::SymbolicSet(mgr, exog_dim,e_lb,e_ub,e_eta);


  tt.tic();
  /* 
  Intermediate Variables representing maxima of small sets of variables
  */
  std::cout << "\nIntermediate Variables" << std::endl;
  std::vector<scots::SymbolicSet> ss_intermed; // intermediate layers ll tree
  ss_intermed.resize(N-2);

  std::vector<scots::FunctionDependency>inter_deps;
  inter_deps.resize(N-1); 

  BDD interconnection = mgr.bddOne();
  tt.tic();
  for (int i = 0; i < N-1; i++){
    std::cout << i << std::endl;
    std::array<double, 1> inter_lb, inter_ub, inter_eta;
    inter_lb[0] = 0;
    inter_ub[0] = 31;
    inter_eta[0] = 1.0;
    if (i < N-2){
      ss_intermed[i] = scots::SymbolicSet(mgr, 1, inter_lb, inter_ub, inter_eta);
    }

    if (i == 0){ // takes two states, outputs to intermed
      inter_deps[i] = scots::FunctionDependency({ss_pre[i], ss_pre[i+1]}, {ss_intermed[i]});
      inter_deps[i].set_dependency(ss_intermed[i][0], {ss_pre[i][0], ss_pre[i+1][0]});
    }
    else if (i == N - 2){ // takes state and intermed, output to exog
      inter_deps[i] = scots::FunctionDependency({ss_intermed[i-1], ss_pre[i+1]}, { ss_exog });
      inter_deps[i].set_dependency(ss_exog[0], {ss_intermed[i-1][0], ss_pre[i+1][0]});
    }
    else{ // takes state and intermed, output to intermed
      inter_deps[i] = scots::FunctionDependency({ss_intermed[i-1], ss_pre[i+1]}, {ss_intermed[i]});
      inter_deps[i].set_dependency(ss_intermed[i][0], {ss_intermed[i-1][0], ss_pre[i+1][0]});
    }

    scots::FunctionAbstracter<std::array<double, 2>, std::array<double, 1> > layer(inter_deps[i], max2);
    interconnection &= layer.compute_abstraction(mgr);

  }
  tt.toc();

  /* Declare dependencies for individual systems */
  std::vector<scots::FunctionDependency> sysdeps(N, scots::FunctionDependency());
  for (int i = 0; i < N; i++){
    sysdeps[i] = scots::FunctionDependency({ss_pre[i], ss_control[i], ss_exog},{ss_post[i]});
    sysdeps[i].set_dependency(ss_post[i][0], {ss_pre[i][0], ss_control[i][0], ss_exog[0]});
  }

  /*Compute system abstractions using dependencies*/
  std::vector<scots::FunctionAbstracter<input_type, state_type> > abs_comp(N, scots::FunctionAbstracter<input_type, state_type>());
  std::vector<BDD> abs_systems(N, mgr.bddOne());
  BDD composed_systems = mgr.bddOne();
  tt.tic();
  for (int i = 0; i < N; i++){
    abs_comp[i] = scots::FunctionAbstracter<input_type, state_type>(sysdeps[i], sys_overapprox);
    std::cout << "System " << i << " abstraction ";
    composed_systems &= abs_comp[i].compute_abstraction(mgr);
    tt.toc();
  }
  tt.toc();

  std::cout << "Composing Systems" << std::endl;
  tt.tic();
  BDD monolithic = composed_systems & interconnection;
  tt.toc();

}

