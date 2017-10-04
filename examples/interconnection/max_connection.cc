/*
 * vehicle.cc
 *
 *  created: Sep 2017
 *   author: Eric Kim
 */

/*
 * information about this example is given ll
 * http://arxiv.org/abs/1313.03715
 * doi: 10.1109/TAC.2016.2593947
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
const int N = 16;

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
/* THIS IS STUPID, BUT I NEED TO DO THIS BECAUSE OF FUNCTIONABSTRACTER'S TEMPLATE ARGUMENTS*/
void max2(std::array<double,2> ll, std::array<double, 2> ur, std::array<double, 1> &o_ll, std::array<double,1> &o_ur){
  int n = ll.size(); 
  for (int i = 0; i < n/2; i++){o_ll[i] = std::max(ll[2*i], ll[2*i+1]); o_ur[i] = std::max(ur[2*i], ur[2*i+1]);}
}
void max4(std::array<double,4> ll, std::array<double, 4> ur, std::array<double, 2> &o_ll, std::array<double,2> &o_ur){
  int n = ll.size(); 
  for (int i = 0; i < n/2; i++){o_ll[i] = std::max(ll[2*i], ll[2*i+1]); o_ur[i] = std::max(ur[2*i], ur[2*i+1]);}
}
void max8(std::array<double,8> ll, std::array<double, 8> ur, std::array<double, 4> &o_ll, std::array<double,4> &o_ur){
  int n = ll.size(); 
  for (int i = 0; i < n/2; i++){o_ll[i] = std::max(ll[2*i], ll[2*i+1]); o_ur[i] = std::max(ur[2*i], ur[2*i+1]);}
}
void max16(std::array<double,16> ll, std::array<double, 16> ur, std::array<double, 8> &o_ll, std::array<double,8> &o_ur){
  int n = ll.size(); 
  for (int i = 0; i < n/2; i++){o_ll[i] = std::max(ll[2*i], ll[2*i+1]); o_ur[i] = std::max(ur[2*i], ur[2*i+1]);}
}
void max32(std::array<double,32> ll, std::array<double, 32> ur, std::array<double, 16> &o_ll, std::array<double,16> &o_ur){
  int n = ll.size(); 
  for (int i = 0; i < n/2; i++){o_ll[i] = std::max(ll[2*i], ll[2*i+1]); o_ur[i] = std::max(ur[2*i], ur[2*i+1]);}
}
void max64(std::array<double,64> ll, std::array<double, 64> ur, std::array<double, 32> &o_ll, std::array<double,32> &o_ur){
  int n = ll.size(); 
  for (int i = 0; i < n/2; i++){o_ll[i] = std::max(ll[2*i], ll[2*i+1]); o_ur[i] = std::max(ur[2*i], ur[2*i+1]);}
}
void max128(std::array<double,128> ll, std::array<double, 128> ur, std::array<double, 64> &o_ll, std::array<double,64> &o_ur){
  int n = ll.size(); 
  for (int i = 0; i < n/2; i++){o_ll[i] = std::max(ll[2*i], ll[2*i+1]); o_ur[i] = std::max(ur[2*i], ur[2*i+1]);}
}

/* Computes AND of BDDs which are a vector of N identical BDDs. Assumes N is a power of 2.*/
// TODO optimize so that the vector is passed by reference
BDD pow2vectorAND(std::vector<BDD> bdds){
  int N = bdds.size(); 
  std::cout << N << std::endl;
  for(int stride = 2; stride <= N; stride *= 2){
    std::cout << "Merging groups of " << stride << std::endl;
    for (int i = 0; i < N; i += stride){
      std::cout << i << "-th group" << std::endl; 
      bdds[i] &= bdds[i+stride/2];
    }
  }
  return bdds[0];
}

int main() {
  /* to measure time */
  TicToc tt;
  /* cudd manager */
  Cudd mgr;
  mgr.AutodynEnable(CUDD_REORDER_SIFT_CONVERGE);
  // mgr.AutodynEnable(CUDD_REORDER_RANDOM_PIVOT);
  //mgr.SetMaxGrowth(2.5);
  mgr.EnableReorderingReporting();
  //mgr.AutodynDisable();

  /* Dynamics for individual subsystem */ 
  auto dynamics = [](const state_type x, const control_type u, const exog_type w) -> state_type {
    state_type post;
    post[0] = std::min(x[0] + u[0], w[0] + 1);
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
    // abs_comp[i] = scots::FunctionAbstracter<input_type, state_type>(sysdeps[i], sys_overapprox);
    // std::cout << "System " << i << " abstraction ";
    // composed_systems &= abs_comp[i].compute_abstraction(mgr);
    // //abs_systems[i] = abs_comp[i].compute_abstraction(mgr);
    // tt.toc();
  }
  tt.toc();

  // std::cout << "Composing smaller systems " << std::endl;
  // BDD interconnected_sys = mgr.bddOne();
  // tt.tic();
  // interconnected_sys = pow2vectorAND(abs_systems);
  // tt.toc();

  tt.tic();

  /* 
  Intermediate Variables representing maxima of small sets of variables
  */
  std::cout << "\nIntermediate Variables" << std::endl;
  int layers = (int) std::log2(N); // WARNING: the code below breaks if N is not a power of 2
  std::vector<scots::SymbolicSet> ss_inter; // intermediate layers ll tree 
  std::vector<scots::FunctionDependency>inter_deps;
  ss_inter.resize(layers-1); inter_deps.resize(layers); 
  int vars_in_layer = N/2;
  for (int i = 0; i < layers; i++){

    std::array<double, N> inter_lb, inter_ub, inter_eta;
    for (int j = 0; j < vars_in_layer; j++){
      inter_lb[j] = 0;
      inter_ub[j] = 31;
      inter_eta[j] = 1.0;
    }
    if (i < layers - 1){
      ss_inter[i] = scots::SymbolicSet(mgr, vars_in_layer, inter_lb,inter_ub,inter_eta);
    }

    /*Declare Symbolic Sets and dependencies*/
    if (i == 0){ // base layer which maps ss_pre -> first layer of 
      inter_deps[i] = scots::FunctionDependency({ss_pre}, {ss_inter[0]});
      for(int j = 0; j < vars_in_layer; j++){
        inter_deps[i].set_dependency(ss_inter[i][j], {ss_pre[2*j][0], ss_pre[2*j+1][0]});
      }
    }
    else if(i == layers-1){ //top layer feeding into exog
      inter_deps[i] = scots::FunctionDependency({ss_inter[i-1]}, {ss_exog});
      inter_deps[i].set_dependency(ss_exog[0], {ss_inter[i-1][0], ss_inter[i-1][1]});
    }
    else{
      inter_deps[i] = scots::FunctionDependency({ss_inter[i-1]}, {ss_inter[i]});
      for(int j = 0; j < vars_in_layer;j++){
        inter_deps[i].set_dependency(ss_inter[i][j], {ss_inter[i-1][2*j], ss_inter[i-1][2*j+1]});
      }
    } 
    //std::cout << inter_deps[i] << std::endl << std::endl;
    vars_in_layer = vars_in_layer / 2;
  }


  /* Abstracter Declarations. This is required b/c templates are all at compile time */
  std::vector<std::vector<BDD> > interconnection_BDDs; interconnection_BDDs.resize(layers);
  std::vector<BDD> layer_BDDs; layer_BDDs.resize(layers);
  BDD interconnection = mgr.bddOne();
  std::cout << "Abstracting Interconnection Layers" << std::endl;
  tt.tic();
  for(int i = layers-1; i >= 0; i--) {
    switch(i){
      case 0:{
        std::cout << "Abstracting " << i << " layers deep" << std::endl;
        scots::FunctionAbstracter<std::array<double, 2>, std::array<double, 1> > layer1(inter_deps[layers-1],max2);
        tt.tic();
        // interconnection_BDDs[0] = layer1.compute_vector_abstraction(mgr);
        // layer_BDDs[i] = pow2vectorAND(layer1.compute_vector_abstraction(mgr));
        interconnection &= pow2vectorAND(layer1.compute_vector_abstraction(mgr));
        // for (int j = 0; j < 1; j++)
          // interconnection &= layer1.compute_abstraction(mgr, j);
        tt.toc();
        break;
      }
      case 1:{
        std::cout << "Abstracting " << i << " layers deep" << std::endl;
        scots::FunctionAbstracter<std::array<double, 4>, std::array<double, 2> > layer2(inter_deps[layers-2],max4);
        // interconnection_BDDs[i] = layer2.compute_vector_abstraction(mgr);
        // layer_BDDs[i] = pow2vectorAND(layer2.compute_vector_abstraction(mgr));
        interconnection &= pow2vectorAND(layer2.compute_vector_abstraction(mgr));
        // for (int j = 0; j < 2; j++)
        //   interconnection &= layer2.compute_abstraction(mgr, j);
        tt.toc();
        break;
      }
      case 2:{
        std::cout << "Abstracting " << i << " layers deep" << std::endl;
        scots::FunctionAbstracter<std::array<double, 8>, std::array<double, 4> > layer3(inter_deps[layers-3],max8);
        // interconnection_BDDs[i] = layer3.compute_vector_abstraction(mgr);
        // layer_BDDs[i] = pow2vectorAND(layer3.compute_vector_abstraction(mgr));
        interconnection &= pow2vectorAND(layer3.compute_vector_abstraction(mgr));
        // for (int j = 0; j < 4; j++)
        //   interconnection &= layer3.compute_abstraction(mgr, j);
        tt.toc();
        break;
      }
      case 3:{ 
        std::cout << "Abstracting " << i << " layers deep" << std::endl;
        scots::FunctionAbstracter<std::array<double, 16>, std::array<double, 8> > layer4(inter_deps[layers-4],max16);
        // interconnection_BDDs[i] = layer4.compute_vector_abstraction(mgr);
        // layer_BDDs[i] = pow2vectorAND(layer4.compute_vector_abstraction(mgr));
        interconnection &= pow2vectorAND(layer4.compute_vector_abstraction(mgr));
        // for (int j = 0; j < 8; j++) 
        //   interconnection &= layer4.compute_abstraction(mgr, j);
        tt.toc();
        break;
      }
      case 4:{
        std::cout << "Abstracting " << i << " layers deep" << std::endl;
        scots::FunctionAbstracter<std::array<double, 32>, std::array<double, 16> > layer5(inter_deps[layers-5],max32);
        // interconnection_BDDs[i] = layer5.compute_vector_abstraction(mgr);
        // layer_BDDs[i] = pow2vectorAND(layer5.compute_vector_abstraction(mgr));
        interconnection &= pow2vectorAND(layer5.compute_vector_abstraction(mgr));
        // for (int j = 0; j < 16; j++)
        //   interconnection &= layer5.compute_abstraction(mgr, j);
        tt.toc();
        break;
      }
      case 5:{
        std::cout << "Abstracting " << i << " layers deep" << std::endl;
        scots::FunctionAbstracter<std::array<double, 64>, std::array<double, 32> > layer6(inter_deps[layers-6],max64);
//        interconnection_BDDs[i] = layer6.compute_vector_abstraction(mgr);
        // layer_BDDs[i] = pow2vectorAND(layer6.compute_vector_abstraction(mgr));
        interconnection &= pow2vectorAND(layer6.compute_vector_abstraction(mgr));
        // for (int j = 0; j < 32; j++)
        //   interconnection &= layer6.compute_abstraction(mgr, j);
        tt.toc();
        break;
      }
    } // end switch

    //layer_BDDs[i] = pow2vectorAND(interconnection_BDDs[i]);
  }
  tt.toc();





  // /** Construct abstraction of interconnected system.
  // First construct system that's a subset of X x W x U x X'
  // Then use this system to get a monolithic system X x U x X'
  // with nondeterminism from W, which is constrained by values ll X
  //  **/



  // std::cout << "Applying interconnection relation" << std::endl;
  // for (int i = 0; i < layers; i++){
  //   std::cout << "Layer " << i << std::endl;
  //   interconnected_sys = interconnected_sys & abs_inter[i];
  // }
  // tt.toc(); tt.tic();
  // std::cout << "Abstracting out internal variables" << std::endl;
  // //interconnected_sys = interconnected_sys.ExistAbstract(exog_product.get_cube(mgr));
  // tt.toc();


}

