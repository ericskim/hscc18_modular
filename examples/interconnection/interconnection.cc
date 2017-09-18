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
const int N = 2;

/*
 * data types for the state space elements and input space
 * elements used in uniform grid
 */
using state_type = std::array<double,state_dim>;
using control_type = std::array<double,control_dim>;
using exog_type  = std::array<double,exog_dim>;
using input_type  = std::array<double,input_dim>;

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



int main() {
  /* to measure time */
  TicToc tt;
  /* cudd manager */
  Cudd mgr;
  mgr.AutodynEnable();
  //mgr.AutodynDisable();

  /* Function */ 
  auto dynamics = [](const state_type x, const control_type u, const exog_type w) -> state_type {
    state_type post;
    post[0] = saturate(x[0] + x[1], 0.0, 40.0);
    post[1] = saturate(x[1] + u[0] + .05*w[0], -1.0, 1.0);
    return post;
  };

  /* Takes an input box and computes an overapproximating box. Both are represented with their 
  */
  auto sys_overapprox = [dynamics](const input_type &i_ll, const input_type &i_ur, state_type& o_ll, state_type& o_ur){
    o_ll = dynamics({i_ll[0], i_ll[1]}, {i_ll[2]}, {i_ll[3]});
    o_ur = dynamics({i_ur[0], i_ur[1]}, {i_ur[2]}, {i_ur[3]});
  };

  /* State spaces */
  std::vector<scots::SymbolicSet> ss_pre; ss_pre.resize(N);
  scots::SymbolicSet pre_product = scots::SymbolicSet();
  std::vector<scots::SymbolicSet> ss_post; ss_post.resize(N);
  scots::SymbolicSet post_product = scots::SymbolicSet();
  state_type s_lb={{0, -1}};
  /* upper bounds of the hyper rectangle */
  state_type s_ub={{40, 1}};
  /* grid node distance diameter */
  state_type s_eta={{.2,.05}};
  for (int i = 0; i < N; i++){
    ss_pre[i] = scots::SymbolicSet(mgr, state_dim,s_lb,s_ub,s_eta);
    ss_post[i] = scots::SymbolicSet(mgr, state_dim,s_lb,s_ub,s_eta);

    pre_product = scots::SymbolicSet(pre_product, ss_pre[i]);
    post_product = scots::SymbolicSet(post_product, ss_post[i]);
  }
  ss_pre[0].print_info(1);
  
  /*Input spaces*/
  std::vector<scots::SymbolicSet> ss_control; ss_control.resize(N);
  scots::SymbolicSet control_product = scots::SymbolicSet();
  /* lower bounds of the hyper rectangle */
  control_type i_lb={{-1}};
  /* upper bounds of the hyper rectangle */
  control_type i_ub={{ 1}};
  /* grid node distance diameter */
  control_type i_eta={{.2}};
  for(int i = 0; i < N; i++){
    ss_control[i] = scots::SymbolicSet(mgr, control_dim,i_lb,i_ub,i_eta);
    control_product = scots::SymbolicSet(control_product, ss_control[i]);
  }

  /* Exogenous spaces */
  std::vector<scots::SymbolicSet> ss_exog; ss_exog.resize(N);
  scots::SymbolicSet exog_product = scots::SymbolicSet();
  exog_type e_lb = {{-40}};
  exog_type e_ub = {{40}};
  exog_type e_eta = {{4.0}};
  for (int i = 0; i < N; i++){
    ss_exog[i] = scots::SymbolicSet(mgr, exog_dim,e_lb,e_ub,e_eta);
    exog_product = scots::SymbolicSet(exog_product, ss_exog[i]);
  }

  /* Declare dependencies for individual systems */
  std::vector<scots::FunctionDependency> sysdeps(N, scots::FunctionDependency());
  for (int i = 0; i < N; i++){
    sysdeps[i] = scots::FunctionDependency({ss_pre[i], ss_control[i], ss_exog[i]},{ss_post[i]});
    sysdeps[i].set_dependency(ss_post[i][0], {ss_pre[i][0], ss_pre[i][1]});
    sysdeps[i].set_dependency(ss_post[i][1], {ss_pre[i][1], ss_control[i][0], ss_exog[i][0]});
  }

  /*Compute system abstractions using dependencies*/
  std::vector<scots::FunctionAbstracter<input_type, state_type> > abs_comp(N, scots::FunctionAbstracter<input_type, state_type>());
  std::vector<BDD> abs_systems(N, mgr.bddOne());
  for (int i = 0; i < N; i++){
    abs_comp[i] = scots::FunctionAbstracter<input_type, state_type>(sysdeps[i], sys_overapprox);
    abs_systems[i] = abs_comp[i].compute_abstraction(mgr);
  }

  /* Declare and abstract interconnection. Specialized to N = 2 */
  std::cout << "Constructing Interconnection Abstraction" << std::endl;
  scots::FunctionDependency inter_dep({ss_pre[0], ss_pre[1]},{ss_exog[0], ss_exog[1]});
  for(int i = 0; i < 2; i++){
    inter_dep.set_dependency(ss_exog[i][0], {ss_pre[0][0], ss_pre[1][0]});
  }
  auto inter_overapprox = [ss_pre](const prod_state_type &ll, const prod_state_type& ur, prod_exog_type &o_ll, prod_exog_type &o_ur){
      std::vector<double> eta = ss_pre[0].get_eta();
      prod_state_type center = {{(ll[0] + ur[0])/2.0,  0, (ll[2] + ur[2])/2.0, 0}};
      for (int i = 0; i < 2; i++){
        o_ll[i] = center[state_dim*i] -.5*(center[0] + center[2]) - (1.0)*(eta[0]);
        o_ur[i] = center[state_dim*i] -.5*(center[0] + center[2]) + (1.0)*(eta[0]);
      }
  };
  scots::FunctionAbstracter<prod_state_type, prod_exog_type> inter_abs_comp(inter_dep, inter_overapprox);
  BDD abs_inter = inter_abs_comp.compute_abstraction(mgr);

  /** Construct abstraction of interconnected system.
  First construct system that's a subset of X x W x U x X'
  Then use this system to get a monolithic system X x U x X'
  with nondeterminism from W, which is constrained by values in X
   **/
  std::cout << "Constructing Monolithic System" << std::endl;
  BDD interconnected_sys = abs_inter;
  for (int i = 0; i < N; i++){
    interconnected_sys &= abs_systems[i];
  }
  std::cout << "Abstracting out internal variables" << std::endl;
  interconnected_sys = interconnected_sys.ExistAbstract(exog_product.get_cube(mgr));

  /** Construct Invariant Set on monlithic space **/
  std::cout << "Constructing Invariant Set" << std::endl;
  BDD inv = mgr.bddZero();
  for (int t = 0; t < 40 ; t++){
    /*Function to declare predicate */
    auto inv_predicate = [t, &ss_pre](const abs_type& idx){
      state_type x;
      ss_pre[0].itox(idx,x);
      /* function returns true if cell associated with x is in invariant set  */
      if ( x[0] <= t + 4.0 && x[0] >= t - 4.0)
        return true;
      return false;
    };

    /* All systems must have x[0] state close to t */
    BDD allin = mgr.bddOne();
    for (int i = 0; i < N; i++){
      allin &= ss_pre[i].ap_to_bdd(mgr,inv_predicate);
    }
    /* Union over all t's*/
    inv |= allin;
  }

  /* Controller synthesis over monolithic system */
  std::cout << "Synthesizing Controller" << std::endl;
  BDD X = mgr.bddZero();
  BDD XX =mgr.bddOne();
  scots::EnfPre enf_pre(mgr,interconnected_sys, pre_product, control_product, post_product);
  // /* the controller */
  BDD C = mgr.bddZero();
  // /* BDD cube for existential abstract inputs */
  const BDD U = control_product.get_cube(mgr);
  std::cout << "\nState Space Size: " << pre_product.get_size(mgr,XX) << std::endl;
  std::cout << "Invariant size: " << pre_product.get_size(mgr,inv) << std::endl << std::endl;
  for(int i=1; XX != X; i++) { 
    X = XX;
    C = enf_pre(X); // (state, input) controlled pre pairs
    XX = C.ExistAbstract(U);
    std::cout << i << "-th predecessor size: " << pre_product.get_size(mgr,XX) << std::endl;
    XX = (XX & (inv));// | forcedout; // 
    std::cout << i << "-th winning domain size: " << pre_product.get_size(mgr,XX) << std::endl;
  }
  
  scots::SymbolicSet controller(pre_product,control_product);
  std::cout << "\nWrite controller to controller.scs \n";
  if(write_to_file(mgr, controller, C,"controller"))
    std::cout << "Done. \n";
}

