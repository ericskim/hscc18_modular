/*
 * iter_consensus.cc
 *
 * Recursive computation of the average value amongst all systems
 *
 *  created: Sep 2017
 *   author: Eric Kim
 */


#include <iostream>
#include <array>
#include <cmath>
#include <list>
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
const int N = 3;

const int inter_dim = 1; 

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

/** @brief Variables IDs that the BDD x depends on **/
std::vector< unsigned int > get_support(const Cudd& mgr, const BDD& x){
  return mgr.SupportIndices({x});
}



int main() {
  /* to measure time */
  TicToc tt;
  /* cudd manager */
  Cudd mgr;
  mgr.AutodynEnable(CUDD_REORDER_SIFT_CONVERGE);
  // mgr.AutodynEnable(CUDD_REORDER_RANDOM_PIVOT);
  //mgr.SetMaxGrowth(2.5);
  //mgr.EnableReorderingReporting();
  //mgr.AutodynDisable();

  /* Dynamics for individual subsystem */ 
  auto dynamics = [](const state_type x, const control_type u, const exog_type w) -> state_type {
    state_type post;
    post[0] = logistic_curve(x[0] + u[0] + .03*w[0], 0, 31);
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
  exog_type e_eta = {{1}};
  ss_exog = scots::SymbolicSet(mgr, exog_dim,e_lb,e_ub,e_eta);
  std::cout << "Exogenous Information" << std::endl;
  ss_exog.print_info(1);

  /* Declare dependencies for individual systems */
  std::vector<scots::FunctionDependency> sysdeps(N, scots::FunctionDependency());
  for (int i = 0; i < N; i++){
    sysdeps[i] = scots::FunctionDependency({ss_pre[i], ss_control[i], ss_exog},{ss_post[i]});
    sysdeps[i].set_dependency(ss_post[i][0], {ss_pre[i][0], ss_control[i][0], ss_exog[0]});
  }
  /* Compute sub-system abstractions using dependencies */
  std::vector<scots::FunctionAbstracter<input_type, state_type> > abs_comp(N, scots::FunctionAbstracter<input_type, state_type>());
  std::vector<BDD> abs_systems(N, mgr.bddOne());
  std::vector<BDD> BDD_components; 
  BDD systems = mgr.bddOne(); 
  BDD interconnected_sys = mgr.bddOne();
  for (int i = 0; i < N; i++){
    abs_comp[i] = scots::FunctionAbstracter<input_type, state_type>(sysdeps[i], sys_overapprox);
    tt.tic();
    std::cout << "System " << i << " abstraction ";
    //abs_systems[i] = abs_comp[i].compute_abstraction(mgr);
    BDD_components.push_back(abs_comp[i].compute_abstraction(mgr));
    tt.toc();
  }


  /* 
  Declare dependencies and abstract interconnection. 
  */
  std::vector<scots::SymbolicSet> ss_intermed; // intermediate layers ll tree
  ss_intermed.resize(N-2);

  std::vector<scots::FunctionDependency>inter_deps; inter_deps.resize(N-1);
  scots::SymbolicSet intermed_product = scots::SymbolicSet();
  BDD abs_inter = mgr.bddOne();
  std::vector<BDD> BDD_recursive_avg(N-1, mgr.bddOne());
  tt.tic();
  for (int i = 0; i < N-1; i++){
    std::cout << i << std::endl;
    std::array<double, 1> inter_lb, inter_ub, inter_eta;
    inter_lb[0] = 0;
    inter_ub[0] = 31;
    inter_eta[0] = 1.0;
    if (i < N-2){
      ss_intermed[i] = scots::SymbolicSet(mgr, 1, inter_lb, inter_ub, inter_eta);
      intermed_product = scots::SymbolicSet(intermed_product, ss_intermed[i]);
    }

    /* 
    Computes recursive average. Index i is captured because it lets the user known how over how many variables the average was computed.
    */
    auto recursive_avg = [i](std::array<double,2> ll, std::array<double, 2> ur, std::array<double, 1> &o_ll, std::array<double,1> &o_ur){
      if (i == 0){
        o_ll[0] = .5*(ll[0] + ll[1]); 
        o_ur[0] = .5*(ur[0] + ur[1]); 
      }
      else{ // left hand side is always the current average of previous i+2 numbers
        o_ll[0] = ((i+1.0)*(ll[0]) + ll[1])/(i+2.0);
        o_ur[0] = ((i+1.0)*(ur[0]) + ur[1])/(i+2.0);
      }
    }; 

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

    scots::FunctionAbstracter<std::array<double, 2>, std::array<double, 1> > layer(inter_deps[i], recursive_avg);
//    BDD_recursive_avg[i] = layer.compute_abstraction(mgr);
    BDD_components.push_back(layer.compute_abstraction(mgr));

  }
  tt.toc();
  std::cout << "Intermediate Product Information" << std::endl; 
  intermed_product.print_info(1);

  /** Construct Invariant Set on monlithic space **/
  std::cout << "Constructing Target Set" << std::endl;
  BDD target = mgr.bddZero();
  for (int t = 0; t < 31 ; t++){
    /* Predicate functions */
    auto inv_predicate = [t, &ss_pre](const abs_type& idx){
      state_type x;
      ss_pre[0].itox(idx,x);
      /* function returns true if cell associated with x is in invariant set */
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

  // /* Handle different types of latent variables */
  BDD X , XX, C = mgr.bddZero(); 
  scots::SymbolicSet latent_product = scots::SymbolicSet(ss_exog, intermed_product); // exogenous and intermediate variables
  scots::SymbolicSet controller(pre_product,control_product);
  const BDD U = control_product.get_cube(mgr);
  const BDD E = latent_product.get_cube(mgr); // all latent variables, exog + intermediate ones
  std::vector<scots::SymbolicSet> ss_latent(ss_intermed);
  ss_latent.emplace(ss_latent.begin(), ss_exog);

  /* Decomposition Predecessor Setup */
  std::cout << "Decomposed Predecessor Setup" << std::endl;
  scots::DecomposedPredecessor decomp_pre(mgr, BDD_components, ss_pre, ss_control, ss_post, ss_latent);
  std::cout << "Decomposed Computing Nonblocking" << std::endl;
  BDD nonblocking = decomp_pre.nonblocking();
  std::cout << "nonblocking state-input: " << (nonblocking) << std::endl;
  print_support(mgr, nonblocking);
  std::cout << "nonblocking states: " << (nonblocking.ExistAbstract(U)) << std::endl;
  print_support(mgr, (nonblocking).ExistAbstract(U));
  // /**
  // Test of is_dependent function. 
  // **/
  // for (size_t i = 0; i<BDD_components.size(); i++){
  //   for (size_t j = 0; j < ss_latent.size(); j++){
  //     std::cout << "Component " << i << " depends on " << j << " " << scots::is_dependent(mgr, BDD_components[i], ss_latent[j].get_cube(mgr)) << std::endl;
  //   }
  // }

  /* Monolithic Predecessor Setup */
  // std::cout << "Monolithic Predecessor Setup" << std::endl;
  // scots::InterconnectedEnfPre enf_pre(mgr,interconnected_sys, pre_product, control_product, post_product, latent_product, abs_inter, abs_systems);
  
  std::cout << std::endl << "State Space Size: " << pre_product.get_size(mgr,mgr.bddOne()) << std::endl;
  std::cout << "Nonblocking Space Size: " << pre_product.get_size(mgr,nonblocking.ExistAbstract(U)) << std::endl;
  std::cout << "Target size: " << pre_product.get_size(mgr,target) << std::endl;
  std::cout << "Nonblocking Pairs: " << controller.get_size(mgr, nonblocking) << std::endl;
  std::cout << "Nonblocking Controls: " << control_product.get_size(mgr, nonblocking) << std::endl  << std::endl;

  /* Controller synthesis over decomposed system */
  BDD dX, dXX, dC = mgr.bddZero();
  tt.tic();
  dX = mgr.bddZero(); dXX = mgr.bddOne();
  for (int i = 0; dXX != dX; i++){
    dX = dXX;
    BDD enf = decomp_pre(dX);
    print_support(mgr, U);
    //std::cout << "enf: " << enf << std::endl;
    print_support(mgr, enf);

    //std::cout << "exists U . enf: " << enf.ExistAbstract(U) << std::endl;
    print_support(mgr, enf.ExistAbstract(U));

    //std::cout << "enf && nonblocking: " << (enf & nonblocking) << std::endl;
    print_support(mgr, (enf & nonblocking));

    //std::cout << pre_product.get_size(mgr, enf.ExistAbstract(U)) << std::endl;
    //std::cout << controller.get_size(mgr, enf) << std::endl;
    dXX = (enf & nonblocking).ExistAbstract(U) & target; // state input pairs
    std::cout << i << "-th winning domain size: " << pre_product.get_size(mgr,dXX) << std::endl << std::endl;
  }
  tt.toc(); 

  // /* Controller synthesis over monolithic system */
  // tt.tic();
  // /*Safety objective*/
  // std::cout<< "Invariance Controller Synthesis" << std::endl;
  // X = mgr.bddZero(); XX = mgr.bddOne();
  // for(int i=1; XX != X; i++) {
  //   X = XX;
  //   C = enf_pre(X); // (state, input) controlled pre pairs
  //   XX = C.ExistAbstract(U*E);
  //   XX = (C & target);
  //   std::cout << i << "-th winning domain size: " << pre_product.get_size(mgr,XX) << std::endl;
  // }

  
  // if(write_to_file(mgr, controller, C.ExistAbstract(E),"consensus_inv_controller"))
  //   std::cout << "Done. \n";

  // auto prod_dynamics = [](prod_state_type &x,  prod_control_type u) {
  //   double avg = 0, w;
  //   for (int i = 0; i < N; i++){
  //     avg += x[i];
  //   }
  //   avg = avg / N;
  //   std::cout << "Average w: " << avg << std::endl;
  //   for (int i = 0; i < N; i++){
  //     w = x[i] - avg;
  //     x[i] = logistic_curve(x[i] + u[i] + .1*w, 0, 31);
  //   }
  // };

  // prod_state_type x={15.4, 15.3, 15, 16, 22, 17};
  // //14.6 15.4 15 16.2 17.1 24.1 
  // int u_index;
  // bool active_control = true;
  // while(true){
  //   std::cout << "Enter an initial 6D state" << std::endl;
  //   std::cin >> x[0] >> x[1] >> x[2] >> x[3] >> x[4] >> x[5];

  //   for(int i=0; i<30; i++) {
  //   //   // returns a std vector with the valid control inputs     
  //     std::cout << "State: ";
  //     for(int j = 0; j < N; j++){
  //       std::cout << x[j] << " ";
  //     }
  //     std::cout<< std::endl;

  //     if (active_control){
  //       std::cout << "Getting Control Input" << std::endl;
  //       auto u = controller.restriction<prod_state_type>(mgr,C,x);
  //       if (u.size() == 0){
  //         std::cout << "No valid control" << std::endl;
  //         break;
  //       }
  //       u_index = u.size() - (rand() % (u.size()/N))*N;//rand() % (u.size()/control_dim);
  //       std::cout << u.size() << std::endl;
  //       std::cout << "Input Index: " << u_index << std::endl;
  //       std::cout << "Input: ";
  //       for(int j = 0; j < N; j++){
  //         std::cout << u[u_index+j] << " ";
  //       }
        
  //       prod_dynamics(x,{u[u_index],u[u_index+1],u[u_index+2],u[u_index+3],u[u_index+4],u[u_index+5]});
  //       std::cout<< std::endl << std::endl;
  //     }
  //     else{
  //       prod_dynamics(x, {0,0,0,0,0,0});
  //     }
  //   }
  // }

}

