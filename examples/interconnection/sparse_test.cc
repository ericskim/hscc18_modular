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
const int state_dim=2;
/* input space dim */
const int input_dim=1;

/*
 * data types for the state space elements and input space
 * elements used in uniform grid and ode solvers
 */
using state_type = std::array<double,state_dim>;
using input_type = std::array<double,input_dim>;

/* abbrev of the type for abstract states and inputs */
using abs_type = scots::abs_type;

/* Identity function */ 
auto identity_post = [](state_type &, const input_type &) {
  /* No changes to state */
};

/* we integrate the growth bound by 0.3 sec (the result is stored in r)  */
auto radius_post = [](state_type &r, const state_type &, const input_type &) {
//  double c = std::abs(u[0])*std::sqrt(std::tan(u[1])*std::tan(u[1])/4.0+1);
  r[0] = .2;
  r[1] = .11;
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
  if(!scots::read_from_file(mgr,ss_pre,"state_pre") ||
     !scots::read_from_file(mgr,ss_post,"state_post")){
   /* lower bounds of the hyper rectangle */
    state_type s_lb={{0,0}};
    /* upper bounds of the hyper rectangle */
    state_type s_ub={{5,10}};
    /* grid node distance diameter */
    state_type s_eta={{.1,.1}};
    /* construct SymbolicSet with the UniformGrid information for the state space
     * and BDD variable IDs for the pre */
    ss_pre = scots::SymbolicSet(mgr, state_dim,s_lb,s_ub,s_eta);
    /* construct SymbolicSet with the UniformGrid information for the state space
     * and BDD variable IDs for the post */
    ss_post = scots::SymbolicSet(mgr, state_dim,s_lb,s_ub,s_eta);
 
    scots::write_to_file(ss_pre,"state_pre");
    scots::write_to_file(ss_post,"state_post");
  }
  std::cout << "Unfiorm grid details:" << std::endl;
  ss_pre.print_info(1);

  /* try to read data from files */
  scots::SymbolicSet ss_input;
  if(!scots::read_from_file(mgr,ss_input,"input_alphabet")) {
    /* construct grid for the input space */
    /* lower bounds of the hyper rectangle */
    input_type i_lb={{-1}};
    /* upper bounds of the hyper rectangle */
    input_type i_ub={{ 1}};
    /* grid node distance diameter */
    input_type i_eta={{.3}};
    ss_input = scots::SymbolicSet(mgr, input_dim,i_lb,i_ub,i_eta);
    scots::write_to_file(ss_input,"input_alphabet");
  }
  ss_input.print_info(1);

  /* Declare dependency and print out */
  scots::DT_Dependency dep(state_dim,input_dim);
  std::vector<std::vector<int>> state_rhs(state_dim, std::vector<int>()), input_rhs(state_dim, std::vector<int>());
  state_rhs[0].push_back({0});
  // input_rhs[0].push_back({0});

  state_rhs[1].push_back({1});
  // input_rhs[1].push_back({0});
  dep.set_rhs(state_rhs, input_rhs);
  std::cout << dep;

  scots::SymbolicSet set;

  /* initialize SymbolicModel class with the abstract state and input alphabet */
  scots::SymbolicModel<state_type,input_type> sym_model(ss_pre,ss_input,ss_post);

  /** Compute transition using sparsity-aware algorithm **/
  BDD sparse_TF;
  if(!scots::read_from_file(mgr,set,sparse_TF,"sparse_tf")) {

    set = scots::SymbolicSet(scots::SymbolicSet(ss_pre,ss_input),ss_post);

    std::cout << "\nComputing the sparse transition function: " << std::endl;
    tt.tic();
    size_t no_trans;
    sparse_TF = sym_model.compute_sparse_gb(mgr,identity_post,radius_post,dep,no_trans);
    tt.toc();
    std::cout << "Number of transitions: " << no_trans << std::endl;
    if(!getrusage(RUSAGE_SELF, &usage))
      std::cout << "Memory per transition: " << usage.ru_maxrss/(double)no_trans << std::endl;
    scots::write_to_file(mgr,set,sparse_TF,"sparse_tf");
  }

  /** Compute transition using brute force algorithm **/
  BDD TF;
  /* does there exist the transition function file ?*/
  if(!scots::read_from_file(mgr,set,TF,"tf")) {
    std::cout << "\nComputing the transition function: " << std::endl;
    tt.tic();
    size_t no_trans;
    TF = sym_model.compute_gb(mgr,identity_post, radius_post,no_trans);
    tt.toc();
    std::cout << "Number of transitions: " << no_trans << std::endl;
    if(!getrusage(RUSAGE_SELF, &usage))
      std::cout << "Memory per transition: " << usage.ru_maxrss/(double)no_trans << std::endl;
    scots::write_to_file(mgr,set,TF,"tf");
  }

  if(sparse_TF == TF){
    std::cout << "\nNo difference between sparse and brute force abstractions!\n"; 
  }
  else {
    std::cout << "\nSparse and brute force abstractions are different!!!\n";
  }

}

