/*
 * test.cc
 *
 *  created on: 20.01.2016
 *      author: rungger
 */

/*
 * a computation to compare the RK4 ode solver with high precision tides solver
 */



#include <iostream>
#include <array>
#include <cmath>

#include "UniformGrid.hh"
#include "AbstractionGB.hh"

#include "TicToc.hh"


extern "C" {
  #include "dp_tides.h"
  #include "aircraftODE.h"
}


/* state space dim */
#define sDIM 3
#define iDIM 2

/* data types for the ode solver */
typedef std::array<double,3> state_type;
typedef std::array<double,2> input_type;

/* we integrate the aircraft ode by 0.25 sec (the result is stored in x)  */
auto  aircraft_post_hi = [](state_type &x, input_type &u) -> void {
	double tolrel = 1.e-16 ;
	double tolabs = 1.e-16;

  int nipt = 1;

	double dt = 0.25; /* (tend-tini)/nipt */

	dp_tides_delta(aircraftODE, NULL, sDIM, iDIM, 0, &x[0], &u[0], 0, dt, nipt, tolrel, tolabs, NULL, NULL);

};

/* forward declaration of the ode solver */
template<class F>
void ode_solver(F rhs, state_type &x, input_type &u, size_t nint, double h);

/* we integrate the vehicle ode by 0.3 sec (the result is stored in x)  */
double mg = 60000.0*9.81;
double mi = 1.0/60000;
auto  aircraft_post = [](state_type &x, input_type &u) -> void {
  /* the ode describing the aircraft */
  auto rhs =[](state_type& xx,  const state_type &x, input_type &u) -> void {
    double c=(1.25+4.2*u[1]);
    xx[0] =        mi*(u[0]*std::cos(u[1])-(2.7+3.08*c*c)*x[0]*x[0]-mg*std::sin(x[1]));
    xx[1] = (1.0/(60000*x[0]))*(u[0]*std::sin(u[1])+68.6*c*x[0]*x[0]-mg*std::cos(x[1]));
    xx[2] = x[0]*std::sin(x[1]);
  };

  size_t nint=5; /* number of intermediate step size */
  double h=0.05; /* h* nint = sampling time */
  ode_solver(rhs,x,u,nint,h);
};



int main() {
  /* to measure time */
  TicToc tt;

  /****************************************************************************/
  /* construct grid for the state space */
  /****************************************************************************/
  /* setup the workspace of the synthesis problem and the uniform grid */
  state_type z={{0.0125,0.0025,0.05}};   
  /* grid node distance diameter */
  state_type eta={{5.0/42,M_PI/12600,4.0/15}};   
  /* lower bounds of the hyper rectangle */
  state_type lb={{58,-3*M_PI/180,0}};  
  /* upper bounds of the hyper rectangle */
  state_type ub={{83,0,56}}; 
  for(size_t i=0; i< sDIM; i++) {
    lb[i]=lb[i]+eta[i]/2+z[i];
    ub[i]=ub[i]-eta[i]/2-z[i];
  }
  scots::UniformGrid<state_type> ss(sDIM,lb,ub,eta,z);
  std::cout << "Unfiorm grid details:" << std::endl;
  //ss.printInfo(1);

  /****************************************************************************/
  /* construct grid for the input space */
  /****************************************************************************/
  /* lower bounds of the hyper rectangle */
  input_type ilb={{0,0}};  
  /* upper bounds of the hyper rectangle */
  input_type iub={{36000,8*M_PI/180}};  
  /* grid node distance diameter */
  input_type ieta={{36000,8.0/9.0*M_PI/180}};  
  scots::UniformGrid<input_type> is(iDIM,ilb,iub,ieta);
  //is.printInfo(1);

 

  size_t N=ss.getN();
  size_t M=is.getN();

   /* state and input variables */
  state_type x;
  state_type y;
  input_type u;

 
  tt.tic();

  /* loop over all states */
  for(size_t i=0; i<N; i+=1e2) {
    for(size_t j=0; j<M; j+=1) {

      /* current state */
      ss.itox(i,x);
      ss.itox(i,y);
      /* current input */
      is.itox(j,u);

      aircraft_post_hi(x,u);
      aircraft_post(y,u);


      for(size_t k=0; k<sDIM; k++) {
        if(std::abs(y[k]-x[k])>1e-8) {
          std::cout << x[0] << " " << x[1] << " " << x[2] << std::endl;
          std::cout << y[0] << " " << y[1] << " " << y[2] << std::endl;
          std::cout << std::endl;
        }

      }

     
 
    }
  }


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


