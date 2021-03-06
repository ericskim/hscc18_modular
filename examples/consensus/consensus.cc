/*
 * consensus.cc
 * 
 * Abstract N bimodal, scalar systems and apply an interconnection where each system 
 * has access to the average state. 
 *
 * Synthesize for a consensus objective. 
 * 
 *  created: Sep 2017
 *  author: Eric Kim
 */

#include <iostream>
#include <algorithm>
#include <array>
#include <cmath>

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
const int N = 6;

const int inter_dim = 2; 

/*
 * data types for the state space elements and input space
 * elements used in uniform grid
 */
using state_type = std::array<double,state_dim>;
using control_type = std::array<double,control_dim>;
using exog_type  = std::array<double,exog_dim>;
using input_type  = std::array<double,input_dim>;
using inter_type = std::array<double, inter_dim>;

/* Data types for the interconnection relation */
using prod_state_type = std::array<double, state_dim * N>;
using prod_control_type = std::array<double, control_dim * N>;
using prod_exog_type = std::array<double, exog_dim * N>;
using prod_intermed_type_lvl2 = std::array<double, inter_dim>;

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

inline double logistic_curve(double x, double lb, double ub, double B = .2){
  double mid = (ub + lb) /2.0;
  double denom = (1 + std::exp(-B*(x-mid)));
  double num = ub - lb;
  return lb + num/denom; 
}

void print_support(const Cudd& mgr, const BDD& x){
  std::vector< unsigned int >  indices = mgr.SupportIndices({x});
  for (size_t i = 0; i < indices.size(); i++){
    std::cout << indices[i] << " ";
  }
  std::cout << std::endl;
}

int main() {
  /* to measure time */
  TicToc tt;
  /* cudd manager and BDD reordering options */
  Cudd mgr;
  mgr.AutodynEnable(CUDD_REORDER_SIFT_CONVERGE);
  // mgr.AutodynEnable(CUDD_REORDER_RANDOM_PIVOT);
  mgr.EnableReorderingReporting(); 

  double K = .1;
  /* Dynamics for individual subsystem */ 
  auto dynamics = [K](const state_type x, const control_type u, const exog_type w) -> state_type {
    state_type post;
    post[0] = logistic_curve(x[0] + u[0] + K*w[0], 0, 31);
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
  control_type i_lb={{-2}};
  /* upper bounds of the hyper rectangle */
  control_type i_ub={{ 2}};
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
  exog_type e_eta = {{2}};
  ss_exog = scots::SymbolicSet(mgr, exog_dim,e_lb,e_ub,e_eta);

  /* Declare dependencies for individual systems */
  std::vector<scots::FunctionDependency> sysdeps(N, scots::FunctionDependency());
  for (int i = 0; i < N; i++){
    sysdeps[i] = scots::FunctionDependency({ss_pre[i], ss_control[i], ss_exog},{ss_post[i]});
    sysdeps[i].set_dependency(ss_post[i][0], {ss_pre[i][0], ss_control[i][0], ss_exog[0]});
  }

  /*Compute system abstractions using dependencies*/
  std::vector<scots::FunctionAbstracter<input_type, state_type> > abs_comp(N, scots::FunctionAbstracter<input_type, state_type>());
  std::vector<BDD> abs_systems(N, mgr.bddOne());
  BDD systems = mgr.bddOne(); 
  for (int i = 0; i < N; i++){
    abs_comp[i] = scots::FunctionAbstracter<input_type, state_type>(sysdeps[i], sys_overapprox);
    tt.tic();
    std::cout << "System " << i << " abstraction ";
    abs_systems[i] = abs_comp[i].compute_abstraction(mgr);
    tt.toc();
  }


  /* 
  Intermediate Variables representing sums of small sets of variables
  */
  std::cout << "\nIntermediate Variables" << std::endl;
  inter_type inter_lb = {{0,0}};
  inter_type inter_ub = {{31,31}};
  inter_type inter_eta = {{2, 2}};
  scots::SymbolicSet ss_inter = scots::SymbolicSet(mgr, inter_dim, inter_lb,inter_ub,inter_eta);
  ss_inter.print_info(1);
  tt.tic();

  /*Interconnection Level 1. Calculate two averages of three variables each. */
  std::cout << "Level 1" << std::endl;
  scots::FunctionDependency intermed_dep_1({pre_product}, {ss_inter});
  intermed_dep_1.set_dependency(ss_inter[0], {ss_pre[0][0], ss_pre[1][0], ss_pre[2][0]});
  intermed_dep_1.set_dependency(ss_inter[1], {ss_pre[3][0], ss_pre[4][0], ss_pre[5][0]});
  std::cout << intermed_dep_1 << std::endl;
  auto inter_overapprox_1 = [ss_pre](const prod_state_type ll, const prod_state_type ur, inter_type &o_ll, inter_type &o_ur){
    std::vector<double> eta = ss_pre[0].get_eta();
    prod_state_type center; 
    // Compute center of box from ll and ur corners 
    for (int i = 0; i < N; i++){
      center[i]   = (ll[i] + ur[i])/2.0;
    }
    o_ll[0] = std::max((1.0/3)*(center[0] + center[1] + center[2] - 1.5*eta[0]),0.0);
    o_ur[0] = std::min((1.0/3)*(center[0] + center[1] + center[2] + 1.5*eta[0]),31.0);
    o_ll[1] = std::max((1.0/3)*(center[3] + center[4] + center[5] - 1.5*eta[0]),0.0); 
    o_ur[1] = std::min((1.0/3)*(center[3] + center[4] + center[5] + 1.5*eta[0]),31.0);
  };
  scots::FunctionAbstracter<prod_state_type, inter_type> inter_abs_1(intermed_dep_1 , inter_overapprox_1);
  BDD abs_inter_1 = inter_abs_1.compute_abstraction(mgr);

  /*Interconnection Level 2. Average of the two averages from Level 1.*/
  std::cout << "\nLevel 2" << std::endl;
  scots::FunctionDependency intermed_dep_2({ss_inter},{ss_exog});
  intermed_dep_2.set_dependency(ss_exog[0], {ss_inter[0], ss_inter[1]});
  std::cout << intermed_dep_2 << std::endl;
  auto inter_overapprox_2 = [pre_product, ss_inter](const prod_intermed_type_lvl2 ll, const prod_intermed_type_lvl2 ur, exog_type &o_ll, exog_type &o_ur){
    std::vector<double> eta = pre_product.get_eta();
    std::vector<double> mu = ss_inter.get_eta();
    exog_type center;
    for (int i = 0; i < 2; i++){
      center[i]  = (ll[i] + ur[i])/2.0;
    }
    o_ll[0] = std::max(0.5*(center[0] + center[1]) - .5*eta[0], 0.0); 
    o_ur[0] = std::min(0.5*(center[0] + center[1]) + .5*eta[0], 31.0);
  };
  scots::FunctionAbstracter<prod_intermed_type_lvl2, exog_type> inter_abs_2(intermed_dep_2 , inter_overapprox_2);
  BDD abs_inter_2 = inter_abs_2.compute_abstraction(mgr); 


  /* 
  Declare and abstract interconnection. 
  */

  BDD abs_inter = abs_inter_1 & abs_inter_2;
  print_support(mgr, abs_inter);
  std::cout << (int)(abs_inter == mgr.bddZero()) << std::endl;
  tt.toc();

  /** Construct abstraction of interconnected system.
  First construct system that's a subset of X x W x U x X'
  Then use this system to get a monolithic system X x U x X'
  with nondeterminism from W, which is constrained by values in X
   **/
  std::cout << "Composing smaller systems " << std::endl;
  BDD interconnected_sys = mgr.bddOne();
  tt.tic();
  for (int i = 0; i < N; i++){
    std::cout << i << std::endl;
    interconnected_sys = interconnected_sys & abs_systems[i];
  }
  tt.toc();

  /** Construct Invariant Set on monlithic space **/
  std::cout << "Constructing Target Set" << std::endl;
  BDD target = mgr.bddZero();
  for (int t = 0; t < 31 ; t++){
    /* Predicate functions */
    auto inv_predicate = [t, &ss_pre](const abs_type& idx){
      state_type x;
      ss_pre[0].itox(idx,x);
      /* function returns true if cell associated with x is in invariant set  */
      if ( x[0] <= t + 3.0 && x[0] >= t - 3.0)
        return true;
      return false;
    };

    /* All systems must have x[0] state close to t */
    BDD allin = mgr.bddOne();
    for (int i = 0; i < N; i++){
      allin &= ss_pre[i].ap_to_bdd(mgr,inv_predicate);
    }
    /* Union over all t's*/
    target |= allin;
  }

  /* Controller synthesis over monolithic system */
  scots::SymbolicSet aux_product = scots::SymbolicSet(ss_exog, ss_inter); // exogenous and intermediate variables
  scots::InterconnectedEnfPre enf_pre(mgr,interconnected_sys, pre_product, control_product, post_product, aux_product, abs_inter, abs_systems);
  scots::SymbolicSet controller(pre_product,control_product);
  // /* the controller */
  BDD X , XX, C = mgr.bddZero(), newbasin;
  // /* BDD cube for existential abstract inputs */
  const BDD U = control_product.get_cube(mgr);
  const BDD E = aux_product.get_cube(mgr); // all auxiliary variables, exog + intermediate ones
  std::cout << "\nState Space Size: " << pre_product.get_size(mgr,mgr.bddOne()) << std::endl;
  std::cout << "Target size: " << pre_product.get_size(mgr,target) << std::endl << std::endl;
  tt.tic();
  /*Safety objective*/
  std::cout<< "Invariance Controller Synthesis" << std::endl;
  X = mgr.bddZero(); XX = mgr.bddOne();
  for(int i=1; XX != X; i++) { 
    X = XX;
    C = enf_pre(X); // (state, input) controlled pre pairs
    XX = C.ExistAbstract(U*E);
    XX = (X & XX & target) & abs_inter.ExistAbstract(E);
    std::cout << i << "-th winning domain size: " << pre_product.get_size(mgr,XX) << std::endl;
  }
  
  if(write_to_file(mgr, controller, C.ExistAbstract(E),"consensus_inv_controller"))
    std::cout << "Done writing controller to file. \n";

  BDD inv = XX;
  std::cout << "Invariant Set BDD Var ID Support" <<std::endl;
  print_support(mgr,inv.ExistAbstract(U*E));
  std::cout << "Target Set BDD Var ID Support" <<std::endl;
  print_support(mgr,target);
  std::cout << "Controller BDD Var ID Support" << std::endl;
  print_support(mgr,C);
  std::cout << "Inter BDD Var ID Support" << std::endl;
  print_support(mgr, abs_inter); 

  /*Reach objective*/
  std::cout<< "Reachability Controller Synthesis" << std::endl;
  X = mgr.bddOne(); XX = mgr.bddZero();
  for(int i = 1; XX != X; i++){
    std::cout << i << "-th reach basin size: " << pre_product.get_size(mgr,XX) << std::endl;
    X = XX;
    XX = (enf_pre(X) | inv) & abs_inter.ExistAbstract(E);
    std::cout << i << "-th predecessor " << pre_product.get_size(mgr,XX) << std::endl;
    newbasin = pre_product.get_grid_bdd(mgr) & XX & (!(C.ExistAbstract(U*E)));
    std::cout << i << "-th new states " << pre_product.get_size(mgr,newbasin) << std::endl << std::endl;
    XX = XX.ExistAbstract(U*E);
    C = C | newbasin;
  }
  tt.toc();

  /* Print final reach set */
  if (false){
    std::ofstream file;
    file.open("better_consensus_reachable.txt");
    auto a = pre_product.bdd_to_grid_points(mgr, XX);
    for(size_t j = 0; j < a.size(); j++){
      file << a[j] << " ";
      if (j % (state_dim *N) == (state_dim *N)-1)
        file << "\n";
    }
    file.close();
  }

  std::cout << "\nWrite controller to better_consensus_controller.scs \n";
  if(write_to_file(mgr, controller, C.ExistAbstract(E),"better_consensus_controller"))
    std::cout << "Done. \n";

  /* Dynamics for monolithic system */
  auto prod_dynamics = [K](prod_state_type &x,  prod_control_type u) {
    double avg = 0, w;
    for (int i = 0; i < N; i++){
      avg += x[i];
    }
    avg = avg / N;
    //std::cout << "Average w: " << avg << std::endl;
    for (int i = 0; i < N; i++){
      w = x[i] - avg;
      x[i] = logistic_curve(x[i] + u[i] + K*w, 0, 31);
    }
  };

  /*
  Simulate Dynamics with synthesized controller 

  Two options:
  - active_control == true: Use synthesized controller
  - active_control == false: All inputs u are equal to zero 
  */
  std::cout << "Simulating Active Controller" << std::endl;
  prod_state_type x={14.6, 15.4, 15, 16.2, 17.1, 24.1};
  int u_index;

  std::ofstream file;
  file.open("traj_active.txt");

  for(int i=0; i<30; i++) {
    file << "State: ";
    for(int j = 0; j < N; j++){
      file << x[j] << " ";
    }
    file<< std::endl;

    file << "Getting Control Input" << std::endl;
    auto u = controller.restriction<prod_state_type>(mgr,C,x);
    if (u.size() == 0){
      file << "No valid control. Exiting" << std::endl;
      break;
    }
    u_index = u.size() - (rand() % (u.size()/N))*N;//rand() % (u.size()/control_dim);
    file << "Number of Permitted Actions: " << u.size() << std::endl;
    file << "Input Index: " << u_index << std::endl;
    file << "Input: ";
    for(int j = 0; j < N; j++){
      file << u[u_index+j] << " ";
    } 
    
    prod_dynamics(x,{u[u_index],u[u_index+1],u[u_index+2],u[u_index+3],u[u_index+4],u[u_index+5]});
    file<< std::endl << std::endl;

  } // close simulation for loop
  file.close();
  std::cout << "Trajectory written to traj_active.txt" << std::endl;

  /* Passive Control */
  std::cout << "Simulating Passive Controller" << std::endl;
  x={14.6, 15.4, 15, 16.2, 17.1, 24.1};
  file.open("traj_passive.txt");
  for(int i=0; i<30; i++) {
    file << "State: ";
    for(int j = 0; j < N; j++){
      file << x[j] << " ";
    }
    file<< std::endl;
    prod_dynamics(x, {0,0,0,0,0,0});
  }
  file.close();
  std::cout << "Trajectory written to traj_passive.txt" << std::endl;



} // close main 

