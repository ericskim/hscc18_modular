/*
 * simulate.cc
 *
 *  created: Jan 2017
 *   author: Matthias Rungger
 */

#include <iostream>
#include <array>
#include <cmath>

/* SCOTS header */
#include "scots.hh"
/* ode solver */
#include "RungeKutta4.hh"

/* state space dim */
const int state_dim=3;
/* input space dim */
const int input_dim=2;

/* sampling time */
const double tau = 0.25;

using state_type = std::array<double,state_dim>;
using input_type = std::array<double,input_dim>;

/* we integrate the aircraft ode by 0.25 sec (the result is stored in x)  */
auto aircraft_post = [] (state_type &x, const input_type &u) {
  /* the ode describing the aircraft */
  auto rhs =[] (state_type& xx,  const state_type &x, const input_type &u) {
    double mg = 60000.0*9.81;
    double mi = 1.0/60000;
    double c=(1.25+4.2*u[1]);
    xx[0] = mi*(u[0]*std::cos(u[1])-(2.7+3.08*c*c)*x[0]*x[0]-mg*std::sin(x[1]));
    xx[1] = (1.0/(60000*x[0]))*(u[0]*std::sin(u[1])+68.6*c*x[0]*x[0]-mg*std::cos(x[1]));
    xx[2] = x[0]*std::sin(x[1]);
  };
  /* use 10 intermediate steps */
  scots::runge_kutta_fixed4(rhs,x,u,state_dim,tau,10);
};

int main() {

  /* define function to check if we are in target */
  auto target = [](const state_type& x) {
    if(         63 <= (x[0]) &&  (x[0]) <=  75 &&
       -3*M_PI/180 <= (x[1]) &&  (x[1]) <=   0 &&
                 0 <= (x[2]) &&  (x[2]) <= 2.5 &&
             -0.91 <= ((x[0]) * std::sin(x[1]) )) {
      return true;
    }
    return false;
  };

  /* read controller from file */
  scots::StaticController con;
  if(!read_from_file(con,"controller")) {
    std::cout << "Could not read controller from controller.scs\n";
    return 0;
  }
  
  std::cout << "\nSimulation:\n " << std::endl;

  state_type x={{81, -1*M_PI/180, 55}};
  while(1) {
    std::vector<input_type> u = con.get_control<state_type,input_type>(x);
    std::cout << x[0] <<  " "  << x[1] << " " << x[2] << "\n";
    //std::cout << u[0][0] <<  " "  << u[0][1] << "\n";
    aircraft_post(x,u[0]);
    if(target(x)) {
      std::cout << "Arrived: " << x[0] <<  " "  << x[1] << " " << x[2] << std::endl;
      break;
    }
  }

  return 1;
}
