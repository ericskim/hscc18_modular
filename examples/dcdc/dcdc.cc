/*
 * Boost Converter 
 *
 *  created on: 28.11.2016
 *      author: rungger
 */

/*
 * information about this example is given in the readme file
 *
 */

#include <iostream>
#include <array>
#include <cmath>

#include "UniformGrid.hh"
#include "AbstractionGB.hh"
#include "SafetyGame.hh"

/* ode solver */
#include "RungeKutta4.hh"

#include "TicToc.hh"
//#include "IO.hh"

/* state space dim */
#define sDIM 2      // dimension of statespace
#define iDIM 1      // dimension of inputspace

/* data types for the ode solver */
typedef std::array<double,sDIM> state_type;
typedef std::array<double,iDIM> input_type;

/* parameters for system dynamics */
double xc=70;
double xl=3;
double rc=0.005;
double rl=0.05;
double ro=1;
double vs=1;

/* parameters for radius calculation */
double k=0.014;
double tau=0.5;			// sampling time
double mu=sqrt(2);

/* ode solver */
OdeSolver ode_solver;


/* we integrate the dcdc ode by 0.5 sec (the result is stored in x)  */
auto system_post = [](state_type &x, const input_type &u) noexcept {
  /* the ode describing the dcdc converter */
  auto rhs =[](state_type& xx,  const state_type &x, const input_type &u) noexcept {
    if(u[0]==1) {
      xx[0]=-rl/xl*x[0]+vs/xl;
      xx[1]=-1/(xc*(ro+rc))*x[1];
    } else {
      xx[0]=-(1/xl)*(rl+ro*rc/(ro+rc))*x[0]-(1/xl)*ro/(5*(ro+rc))*x[1]+vs/xl;
      xx[1]=(1/xc)*5*ro/(ro+rc)*x[0]-(1/xc)*(1/(ro+rc))*x[1];
    }
	};
  ode_solver(rhs,x,u,sDIM,tau);
};

/* we integrate the growth bound by 0.5 sec (the result is stored in r)  */
auto radius_post = [](state_type &r, const state_type&, const input_type &u) noexcept {
  /* the ode for the growth bound */
  auto rhs =[](state_type& rr,  const state_type &r, const input_type &u) noexcept {
    if(u[0]==1) {
      rr[0]=-rl/xl*r[0];
      rr[1]=-1/(xc*(ro+rc))*r[1];
    } else {
      rr[0]=-(1/xl)*(rl+ro*rc/(ro+rc))*r[0]+(1/xl)*ro/(5*(ro+rc))*r[1];
      rr[1]=5*(1/xc)*ro/(ro+rc)*r[0]-(1/xc)*(1/(ro+rc))*r[1];
    }
	};
  ode_solver(rhs,r,u,sDIM,tau);
};

int main() {
  /* to measure time */
  TicToc tt;

  /****************************************************************************/
  /* construct grid for the state space */
  /****************************************************************************/
  /* setup the workspace of the synthesis problem and the uniform grid */
  /* lower bounds of the hyper rectangle */
  //state_type lb={{1.15,5.45}};
  ///* upper bounds of the hyper rectangle */
  //state_type ub={{1.55,5.85}};
  ///* grid node distance diameter */
  //state_type eta={{2/4e3,2/4e3}};
  state_type lb={{1.15,5.45}};
  /* upper bounds of the hyper rectangle */
  state_type ub={{1.55,5.85}};
  /* grid node distance diameter */
  state_type eta={{2/4e3,2/4e3}};
  scots::UniformGrid<state_type> ss(sDIM,lb,ub,eta);
  std::cout << "Unfiorm grid details:" << std::endl;
  ss.printInfo(1);

  /****************************************************************************/
  /* construct grid for the input space */
  /****************************************************************************/
  /* lower bounds of the hyper rectangle */
  input_type ilb={{1}};
  /* upper bounds of the hyper rectangle */
  input_type iub={{2}};
  /* grid node distance diameter */
  input_type ieta={{1}};
  scots::UniformGrid<input_type> is(iDIM,ilb,iub,ieta);

  std::cout << "Compute absraction "<< std::endl;
  /* transition system to be computed */
  scots::TransitionSystem ts;
  tt.tic();
  scots::AbstractionGB<state_type,input_type> abs(ss,is,ts);
  abs.computeTransitionRelation(system_post,radius_post);

  std::cout << "Number of transitions: " << ts.getNoTransitions() << std::endl;
  tt.toc();


  /* calculate maximal fixed point */
  /* define function to check if the cell is in the safe set?  */
  state_type x;
  auto safeset = [&](const size_t idx) noexcept {
    ss.itox(idx,x);
    /* function returns 1 if cell associated with x is in target set  */
    if (lb[0] <= (x[0]-eta[0]/2.0) && (x[0]+eta[0]/2.0)<= ub[0] && lb[1] <= (x[1]-eta[1]/2.0) &&  (x[1]+eta[1]/2.0) <= ub[1])
      return true;
    return false;
  };

  std::cout << "Start synthesis "<< std::endl;
  /* save the result of safety controller in safe.scs */
  scots::SafetyGame safety(ts);
  tt.tic();
  safety.solve(safeset);
  tt.toc();

  std::cout << "Size: " << safety.size() << std::endl;
  std::cout << "Size pairs: " << safety.sizePairs() << std::endl;

  return 1;
}

