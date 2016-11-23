/*
 * unicycle.cc
 *
 *  created on: 26.10.2015
 *      author: rungger
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

/* we integrate the unicycle ode by 0.3 sec (the result is stored in x)  */
auto  unicycle_post = [](state_type &x, input_type &u) -> void {

  /* the ode describing the unicycle */
  auto rhs =[](state_type& xx,  const state_type &x, input_type &u) -> void {
      xx[0] = u[0]*std::cos(x[2]);
      xx[1] = u[0]*std::sin(x[2]);
      xx[2] = u[1];
  };

  size_t nint=5; /* number of intermediate step size */
  double h=0.06; /* h* nint = sampling time */
  ode_solver(rhs,x,u,nint,h);
};

/* computation of the growth bound (the result is stored in r)  */
auto radius_post = [](state_type &r, input_type &u) -> void {
    r[0] = r[0]+r[2]*std::abs(u[0])*0.3;
    r[1] = r[1]+r[2]*std::abs(u[0])*0.3;
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
  state_type ub={{24.4,23.6,M_PI+0.4}};
  /* grid node distance diameter */
  state_type eta={{.1,.1,.1}};
  scots::UniformGrid<state_type> ss(sDIM,lb,ub,eta);
  std::cout << "Unfiorm grid details:" << std::endl;
  ss.printInfo(1);

  /****************************************************************************/
  /* construct grid for the input space */
  /****************************************************************************/
  /* lower bounds of the hyper rectangle */
    input_type ilb={{-0.6,-0.5}};
    /* upper bounds of the hyper rectangle */
    input_type iub={{0.6,0.5}};
    /* grid node distance diameter */
    input_type ieta={{.1,.05}};
  scots::UniformGrid<input_type> is(iDIM,ilb,iub,ieta);

  /****************************************************************************/
  /* set up constraint functions with obtacles */
  /****************************************************************************/
  double H[21][4] = {
    
      {0.0098,4.2688,0,23.72000},
      {6.7181,18.4078,7.0000,15.7800},
      {4.0458,17.9950,18.1600,23.7400},
      {7.0428,24.2,0.00,4.1800},
      {20.8024,24.4226,3.0200,23.5400},
      
      {1.9070-0.3,4.2301+0.3,10.3000-0.3,11.3600+0.3},
      {4.9627,5.8918+0.3,18.3400-0.3,21.9400+0.3},
      {4.6627,5.8918,11.4800-0.3,14.9400+0.3},
      {4.6627,5.8918,6.7400-0.3,9.2800+0.3},
      {4.6627,5.8918,0,3.3200+0.3},
      {11.9323-0.3,12.9195+0.3,7.4200-0.3,14.5800+0.3},
      {11.9323-0.3,12.9195+0.3,18.4200-0.3,22.1400+0.3},
      {11.9323-0.3,12.9195+0.3,0.3000-0.3,3.2200+0.3},
      {7.7150-0.3,10.9093+0.3,4.6400-0.3,5.6000+0.3},
      {14.1426-0.3,17.4336+0.3,4.6400-0.3,5.6000+0.3},
      {7.3085-0.3,10.8706+0.3,16.1600-0.3,17.1000+0.3},
      {14.1232-0.3,17.9757+0.3,16.1600-0.3,17.1000+0.3},
      {18.8245-0.3,19.7536+0.3,7.1600-0.3,10.0200+0.3},
      {18.8245-0.3,19.7536+0.3,12.3600-0.3,15.3600+0.3},
      {18.8245-0.3,19.7536+0.3,18.2600-0.3,21.7400+0.3},
      {20.8992-0.3,23.0674+0.3,10.3800-0.3,11.1600+0.3}

      
      
      
      
  };
    
    /*double H[21][4] = {
        {0+0.2,4.3139-0.2,16.7000+0.2,23.5000-0.2},
        {1.9587+0.2,4.3139-0.2,12.1200+0.2,17.7200-0.2},
        {0.2,1.2938-0.2,12.1200+0.2,17.7200-0.2},
        {6.8181+0.2,11.3997-0.2,13.4600+0.2,15.3000-0.2},
        {6.8374+0.2,8.8442-0.2,10.1400+0.2,14.7600-0.2},
        {6.8568+0.2,11.4190-0.2,6.5800+0.2,11.1400-0.2},
        {9.7221+0.2,11.3029-0.2,19.3600+0.2,23.5000-0.2},
        {6.8762+0.2,11.3610-0.2,18.0400+0.2,20.6000-0.2},
        {13.6134+0.2,24.2-0.2,0.2+0.2,4-0.2},
        {13.8845+0.2,18.2144-0.2,6.5600+0.2,8.6600-0.2},
        {13.8651+0.2,18.1757-0.2,11.1000+0.2,15.2800-0.2},
        {6.8181+0.2,11.3416-0.2,0.2+0.2,3.8600-0.2},
        {21.1798-0.5+0.2,24.2-0.2,12.4600-0.5+0.2,19.4400+0.5-0.2},
        {20.6798+0.2,24.2-0.2,11.9600+0.2,19.9400-0.2},
        {16.5562+0.2,18.1950-0.2,21.1400+0.2,23.5000-0.2},
        {20.6411+0.2,24.2-0.2,6.4200+0.2,9.5000-0.2},
        {13.8070+0.2,18.1563-0.2,17.8600+0.2,22.4000-0.2},
        {14.3070-0.5+0.2,17.6563+0.5-0.2,18.3600-0.5+0.2,21.9000+0.5-0.2},
        {0.2,1.4874-0.2,4.1800+0.2,9.6000-0.2},
        {1.9394+0.2,4.3914-0.2,4.2000+0.2,9.5600-0.2},
        {0.2,4.4301-0.2,2.0800+0.2,5.5000-0.2}
        
    };*/

    
    
    
    

  /* overflow function returns 1 if x \in overflow symbol  */
  auto overflow = [&](const state_type &x) -> bool {
    double c1= eta[0]/2.0+1e-16;
    double c2= eta[1]/2.0+1e-16;
    for(size_t i=0; i<21; i++) {
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

  abs.computeTransitionRelation(unicycle_post, radius_post, overflow);

  tt.toc(); 
  std::cout << "Number of transitions: " << ts.getNoTransitions() << std::endl;

  /* define function to check if the cell is in the  target set?  */
  state_type x;
  auto target = [&](const size_t idx) -> bool {
    ss.itox(idx,x);
    /* function returns 1 if cell associated with x is in target set  */
    if (19 <= (x[0]-eta[0]/2.0) && (x[0]+eta[0]/2.0)<= 20 && 16.5 <= (x[1]-eta[1]/2.0) &&  (x[1]+eta[1]/2.0) <= 17.5)
      return true;
    //if (5 <= (x[0]-eta[0]/2.0) && (x[0]+eta[0]/2.0)<= 5.5 && 5 <= (x[1]-eta[1]/2.0) &&  (x[1]+eta[1]/2.0) <= 5.5)
      //   return true;
  
      
    else
      return false;
  };


  ss.clearAbstractSet();
  ss.addIndices(target);
  scots::IO::writeToFile(&ss,"target.scs");

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


