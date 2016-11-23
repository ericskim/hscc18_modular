/*
 * vehicle_alt.cc
 *
 *  created on: 26.10.2015
 *      author: rungger
 */

/*
 * information about this example is given in the readme file
 *
 */

#include <iostream>
#include <array>

#include "UniformGrid.hh"
#include "TransitionSystem.hh"
#include "AbstractionGB.hh"
#include "ReachabilityGame.hh"

#include "TicToc.hh"
#include "IO.hh"
/* state space dim */
#define sDIM 3
#define iDIM 2

/* data types for the state space elements and input space
 * elements used in uniform grid and ode solvers 
 */
typedef std::array<double,3> state_type;
typedef std::array<double,2> input_type;

/* forward declaration of the ode solver */
template<class F>
void ode_solver(F rhs, state_type &x, input_type &u, size_t nint, double h);

/* we integrate the vehicle ode by 0.3 sec (the result is stored in x)  */
auto  vehicle_post = [](state_type &x, input_type &u) -> void {

  /* the ode describing the vehicle */
  auto rhs =[](state_type& xx,  const state_type &x, input_type &u) -> void {
      double alpha=std::atan(std::tan(u[1])/2.0);
      xx[0] = u[0]*std::cos(alpha+x[2])/std::cos(alpha);
      xx[1] = u[0]*std::sin(alpha+x[2])/std::cos(alpha);
      xx[2] = u[0]*std::tan(u[1]);
  };

  size_t nint=5; /* number of intermediate step size */
  double h=0.06; /* h* nint = sampling time */
  ode_solver(rhs,x,u,nint,h);
};

/* we integrate the growth bound by 0.3 sec (the result is stored in r)  */
auto radius_post = [](state_type &r, input_type &u) -> void {
    double c = std::abs(u[0]*std::sqrt(std::tan(u[1])*std::tan(u[1])/4.0+1));
    r[0] = r[0]+c*r[2]*0.3;
    r[1] = r[1]+c*r[2]*0.3;
};

int main() {
  /* to measure time */
  TicToc tt;

  /****************************************************************************/
  /* construct grid for the state space */
  /****************************************************************************/
  /* setup the workspace of the synthesis problem and the uniform grid */
  /* lower bounds of the hyper rectangle */
  state_type lb={{0,0,-M_PI-0.4}};  
  /* upper bounds of the hyper rectangle */
  state_type ub={{10,10,M_PI+0.4}}; 
  /* grid node distance diameter */
  state_type eta={{.1,.1,.1}};   
  scots::UniformGrid<state_type> ss(sDIM,lb,ub,eta);
  std::cout << "Unfiorm grid details:" << std::endl;
  ss.printInfo(1);

  /****************************************************************************/
  /* construct grid for the input space */
  /****************************************************************************/
  /* lower bounds of the hyper rectangle */
  input_type ilb={{-1,-1}};  
  /* upper bounds of the hyper rectangle */
  input_type iub={{1,1}};  
  /* grid node distance diameter */
  input_type ieta={{.3,.3}};  
  scots::UniformGrid<input_type> is(iDIM,ilb,iub,ieta);

  /****************************************************************************/
  /* set up constraint functions with obtacles */
  /****************************************************************************/
  double H[15][4] = {
    { 1  , 1.2, 0  ,   9 },
    { 2.2, 2.4, 0  ,   5 },
    { 2.2, 2.4, 6  ,  10 },
    { 3.4, 3.6, 0  ,   9 },
    { 4.6, 4.8, 1  ,  10 },
    { 5.8, 6  , 0  ,   6 },
    { 5.8, 6  , 7  ,  10 },
    { 7  , 7.2, 1  ,  10 },
    { 8.2, 8.4, 0  ,  8.5},
    { 8.4, 9.3, 8.3,  8.5},
    { 9.3, 10 , 7.1,  7.3},
    { 8.4, 9.3, 5.9,  6.1},
    { 9.3, 10 , 4.7,  4.9},
    { 8.4, 9.3, 3.5,  3.7},
    { 9.3, 10 , 2.3,  2.5}
  };

  /* overflow function returns 1 if x \in overflow symbol  */
  auto overflow = [&](const state_type &x) -> bool {
    double c1= eta[0]/2.0+1e-16;
    double c2= eta[1]/2.0+1e-16;
    for(size_t i=0; i<15; i++) {
      if ((H[i][0]-c1) <= x[0] && x[0]<= (H[i][1]+c1) && (H[i][2]-c2) <= x[1] && x[1] <= (H[i][3]+c2))
        return true;
    }
    return false;
  };

  ss.addGridPoints(overflow);

  scots::IO::writeToFile(&ss,"obstacles.scs");

  /* transition system to be computed */
  scots::TransitionSystem ts;

  tt.tic();

  scots::AbstractionGB<state_type,input_type> abs(&ss,&is,&ts);

  abs.computeTransitionRelation(vehicle_post, radius_post, overflow);

  tt.toc(); 
  std::cout << "Number of transitions: " << ts.getNoTransitions() << std::endl;

  /* define function to check if the cell is in the  target set?  */
  state_type x;
  auto target = [&](const size_t idx) -> bool {
    ss.itox(idx,x);
    /* function returns 1 if cell associated with x is in target set  */
    if (9 <= (x[0]-eta[0]/2.0) && (x[0]+eta[0]/2.0)<= 9.5 && 0 <= (x[1]-eta[1]/2.0) &&  (x[1]+eta[1]/2.0) <= 0.5) 
      return true;
    else
      return false;
  };

  ss.fillAbstractSet();
  ss.remIndices(target);
  ss.remGridPoints(overflow);
  scots::IO::writeToFile(&ss,"problemdomain.scs");

  tt.tic();

  scots::ReachabilityGame reach(&ts);

  reach.solve(target);

  scots::IO::writeControllerToFile(&reach,"reach.scs",&ss,&is);

  std::cout << "Size: " << reach.sizeOfDomain() << std::endl;
  tt.toc();

  return 1;
}

template<class F>
void ode_solver(F rhs, state_type &x, input_type &u, size_t nint, double h) {
  /* runge kutte order 4 */
  state_type k[4];
  state_type tmp;

  for(size_t t=0; t<nint; t++) {
    rhs(k[0],x, u);
    for(size_t i=0;i<sDIM;i++)
      tmp[i]=x[i]+h/2*k[0][i];

    rhs(k[1],tmp, u);
    for(size_t i=0;i<sDIM;i++)
      tmp[i]=x[i]+h/2*k[1][i];

    rhs(k[2],tmp, u);
    for(size_t i=0;i<sDIM;i++)
      tmp[i]=x[i]+h*k[2][i];

    rhs(k[3],tmp, u);
    for(size_t i=0; i<sDIM; i++)
      x[i] = x[i] + (h/6)*(k[0][i] + 2*k[1][i] + 2*k[2][i] + k[3][i]);
  }
}


