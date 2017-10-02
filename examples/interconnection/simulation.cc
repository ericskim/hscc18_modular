/*
 * simulate.cc
 *
 *  created: Sep 2017
 *   author: Eric Kim 
 */

#include <iostream>
#include <array>
 #include <stdlib.h>

/* SCOTS header */
#include "scots.hh"
/* ode solver */
#include "RungeKutta4.hh"

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
using prod_control_type = std::array<double, control_dim * N>;

/* abbrev of the type for abstract states and inputs */
using abs_type = scots::abs_type;

/* data type used by the ODE solver */
using state_type = std::array<double,state_dim>;

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

/* Monolithic dynamics */ 
auto dynamics = [](prod_state_type& x, const prod_control_type u) {
  double w0 = x[0] - .5*(x[0] + x[2]);
  double w1 = x[2] - .5*(x[0] + x[2]);
  double w2 = x[4] - .5*(x[0] + x[2]);
  x[0] = saturate(x[0] + x[1], 0.0, 32.0);
  x[1] = saturate(x[1] + u[0] + .03*w0, -1.4, 1.4);
  x[2] = saturate(x[2] + x[3], 0.0, 32.0);
  x[3] = saturate(x[3] + u[1] + .03*w1, -1.4, 1.4);
  // x[4] = saturate(x[4] + x[5], 0.0, 40.0);
  // x[5] = saturate(x[5] + u[2] + .03*w2, -1.4, 1.4); 
};

int main() {

  /* Cudd manager */
  Cudd manager;

  /* read controller from file */
  BDD C;
  scots::SymbolicSet con;
  if(!read_from_file(manager,con,C,"controller")) {
    std::cout << "Could not read controller from controller.scs\n";
    return 0;
  }
  std::cout << "Controller Information" << std::endl; 
  con.print_info(1);
  print_support(manager, C);
  std::cout << "\nSimulation:\n " << std::endl;

  prod_state_type x={24.0, 0.4, 22.0, 1.0};//, 12, 0.0};

  for(int i=0; i<50; i++) {
  //   // returns a std vector with the valid control inputs 
     auto u = con.restriction<prod_state_type>(manager,C,x);
     int u_index = u.size() - 10 ;//rand() % (u.size()/2);
     std::cout << u.size() << std::endl;
     std::cout << "State: " << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << "\n";
     if (u.size() == 0){
      std::cout << "No valid control" << std::endl;
      break;
     }
     std::cout << "Input: " << u[u_index] << " " << u[u_index+1] << "\n\n";
     dynamics(x,{u[u_index],u[u_index+1]});
  }

  return 1;
}
