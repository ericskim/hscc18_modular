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
const int state_dim=6;
/* input space dim */
const int control_dim=6;
/* exog space dim*/
const int exog_dim = 1;
/* input space is a cartesian product*/
const int input_dim = state_dim + control_dim + exog_dim; 

/*
 * data types for the state space elements and input space
 * elements used in uniform grid
 */
using state_type = std::array<double,state_dim>;
using control_type = std::array<double,control_dim>;
using exog_type  = std::array<double,exog_dim>;
using input_type  = std::array<double,input_dim>;

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

// 31 / (1 + e^(-.2*(x-16)))
inline double logistic_curve(double x, double lb, double ub, double B = .2){
  double mid = (ub + lb) /2.0;
  double denom = (1 + std::exp(-B*(x-mid)));
  double num = ub - lb;
  return lb + num/denom; 
}

auto dynamics = [](state_type &x, control_type u) {
  double avg = 0, w;
  for (int i = 0; i < state_dim; i++){
    avg += x[i];
  }
  avg = avg / state_dim;
  for (int i = 0; i < state_dim; i++){
    w = x[i] - avg;
    x[i] = logistic_curve(x[i] + u[i] + .1*w, 0, 31);
  }
};

int main() {

  /* Cudd manager */
  Cudd manager;

  /* read controller from file */
  BDD C;
  scots::SymbolicSet con, ss_pre, ss_control;
  bool active_control = false;
  if(active_control && !read_from_file(manager,con,C,"consensus_controller")) {
    std::cout << "Could not read controller from consensus_controller.scs\n";
    return 0;
  }
  if(active_control){
    std::cout << "Controller Information" << std::endl; 
    con.print_info(1);
    ss_pre = scots::SymbolicSet(con, {0,1,2,3,4,5});
    ss_control = scots::SymbolicSet(con, {6,7,8,9,10,11});
    // std::ofstream file;
    // file.open("consensus_reachable.txt");
    // auto a = ss_pre.bdd_to_grid_points(manager, C.ExistAbstract(ss_control.get_cube(manager)));
    // for(size_t j = 0; j < a.size(); j++){
    //   // if (j % (state_dim *N) == 0){
    //   //   file << i << " ";
    //   // }
    //   file << a[j] << " ";
    //   if (j % (state_dim) == (state_dim)-1)
    //     file << "\n";
    // }
    // file.close(); 
    print_support(manager, C);
  }



  //2 12 8 14 8 10
  state_type x={15.5, 15.8, 15, 16, 24, 16};
  int u_index;
  std::cout << "\nSimulation:\n " << std::endl;

      // std::cout << "Enter an initial 6D state" << std::endl;
      // std::cin >> x[0] >> x[1] >> x[2] >> x[3] >> x[4] >> x[5];

      for(int i=0; i<50; i++) {
      //   // returns a std vector with the valid control inputs     
        std::cout << "State: ";
        for(int j = 0; j < state_dim; j++){
          std::cout << x[j] << " ";
        }
        std::cout<< std::endl;

        if (active_control){
          std::cout << "Getting Control Input" << std::endl;
          auto u = con.restriction<state_type>(manager,C,x);
          u_index = u.size() - (rand() % (u.size()/control_dim))*control_dim;//rand() % (u.size()/control_dim);
          std::cout << u.size() << std::endl;
          if (u.size() == 0){
            std::cout << "No valid control" << std::endl;
            break;
          }
          std::cout << "Input Index: " << u_index << std::endl;
          std::cout << "Input: ";
          for(int j = 0; j < control_dim; j++){
            std::cout << u[u_index+j] << " ";
          }
          std::cout<< std::endl << std::endl;

          dynamics(x,{u[u_index],u[u_index+1],u[u_index+2],u[u_index+3],u[u_index+4],u[u_index+5]});
        }
        else{
          dynamics(x, {0,0,0,0,0,0});
        }
      }
    
  



  return 1;
}
