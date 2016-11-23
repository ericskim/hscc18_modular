/*
 * Boost Converter 
 *
 *  created on: 17.12.2015
 *      author: tong
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
#include "TicToc.hh"
#include "SafetyGame.hh"
#include "IO.hh"

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
double tau_s=0.5;			// sampling time
double tau_d=2;
double mu=sqrt(2);

/* we integrate the vehicle ode by 0.3 sec (the result is stored in x)  */
auto  system_post = [](state_type &x, input_type &u) -> void {

  /* the ode describing the vehicle */
  auto rhs =[](state_type& xx,  const state_type &x, input_type &u) -> void {

    switch((int)u[0])
    {
      case 1:
        xx[0]=-rl/xl*x[0]+vs/xl;
        xx[1]=-1/(xc*(ro+rc))*x[1];
        break;
      case 2:
        xx[0]=-(1/xl)*(rl+ro*rc/(ro+rc))*x[0]-(1/xl)*ro/(5*(ro+rc))*x[1]+vs/xl;
        xx[1]=(1/xc)*5*ro/(ro+rc)*x[0]-(1/xc)*(1/(ro+rc))*x[1];
        break;
    }

	};

	/* runge kutte order 4 */
  state_type k[4];
  state_type tmp;

  size_t nint=5; /* number of intermediate step size */
  double h=tau_s/nint; /* h* nint = sampling time */

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
		
};

auto overflow = [](const state_type &) -> bool {
  return false;
};


/* we integrate the growth bound by 0.3 sec (the result is stored in r)  */
auto radius_post = [](state_type &r, input_type &u) -> void {

  /* the ode describing the vehicle */

  auto rhs =[](state_type& rr,  const state_type &r, input_type &u) -> void {

    switch((int)u[0])
		{
      case 1:
        rr[0]=-rl/xl*r[0];
        rr[1]=-1/(xc*(ro+rc))*r[1];
        break;
      case 2:
        rr[0]=-(1/xl)*(rl+ro*rc/(ro+rc))*r[0]+(1/xl)*ro/(5*(ro+rc))*r[1];
        rr[1]=5*(1/xc)*ro/(ro+rc)*r[0]-(1/xc)*(1/(ro+rc))*r[1];
        break;
      }

	};

	/* runge kutte order 4 */
  state_type k[4];
  state_type tmp;

  size_t nint=5; /* number of intermediate step size */
  double h=tau_s/nint; /* h* nint = sampling time */

	for(size_t t=0; t<nint; t++) {
		rhs(k[0],r, u);
		for(size_t i=0;i<sDIM;i++)
		  tmp[i]=r[i]+h/2*k[0][i];

		rhs(k[1],tmp, u);
		for(size_t i=0;i<sDIM;i++)
		  tmp[i]=r[i]+h/2*k[1][i];

		rhs(k[2],tmp, u);
		for(size_t i=0;i<sDIM;i++)
		  tmp[i]=r[i]+h*k[2][i];

		rhs(k[3],tmp, u);
		for(size_t i=0; i<sDIM; i++)
		  r[i] = r[i] + (h/6)*(k[0][i] + 2*k[1][i] + 2*k[2][i] + k[3][i]);
	}

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
//  state_type eta={{1/(4000*sqrt(2)),1/(4000*sqrt(2))}};
  state_type eta={{1/(4000*sqrt(2)),1/(4000*sqrt(2))}};
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

  scots::TransitionSystem ts;
  tt.tic();
  scots::AbstractionGB<state_type,input_type> abs(&ss,&is,&ts);
  abs.computeTransitionRelation(system_post,radius_post,overflow);

  std::cout << "Number of transitions: " << ts.getNoTransitions() << std::endl;

  tt.toc();


  /* calculate maximal fixed point */
  state_type lb_s={{1.15,5.45}};
  state_type ub_s={{1.55,5.85}};

  /* define function to check if the cell is in the safe set?  */
    state_type x;
    auto specification = [&](const size_t idx) -> bool {
      ss.itox(idx,x);
      /* function returns 1 if cell associated with x is in target set  */
      if (lb_s[0] <= (x[0]-eta[0]/2.0) && (x[0]+eta[0]/2.0)<= ub_s[0] && lb_s[1] <= (x[1]-eta[1]/2.0) &&  (x[1]+eta[1]/2.0) <= ub_s[1])
        return true;
      else
        return false;
    };

  /* save the result of safety controller in safe.scs */
  scots::SafetyGame sg(&ts);
  tt.tic();
  sg.solve(specification);
  tt.toc();
  scots::IO::writeControllerToFile(&sg,"safe.scs",&ss,&is);
  std::cout << "Size: " << sg.sizeOfDomain() << std::endl;

  /* save the safety domain in safedomain.scs */
  auto result = [&](const size_t idx) -> bool {
    return sg.ifInDomain(idx);
  };
  ss.addIndices(result);
  scots::IO::writeToFile(&ss,"safedomain.scs");


  return 1;
}

