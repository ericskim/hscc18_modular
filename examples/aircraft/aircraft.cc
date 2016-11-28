/*
 * aircraft.cc
 *
 *  created on: 18.01.2016
 *      author: rungger
 */

/*
 * information about this example is given in
 * http://arxiv.org/abs/1503.03715
 *
 */

#include <iostream>
#include <array>

#include "UniformGrid.hh"
#include "TransitionSystem.hh"
#include "AbstractionGB.hh"
#include "ReachabilityGame.hh"

/* ode solver */
#include "RungeKutta4.hh"

#include "TicToc.hh"

/* state space dim */
#define sDIM 3
#define iDIM 2

/* data types for the state space elements and input space
 * elements used in uniform grid and ode solvers
 */
typedef std::array<double,3> state_type;
typedef std::array<double,2> input_type;

/* sampling time */
const double tau = 0.25;

/* ode solver */
OdeSolver ode_solver;

/* we integrate the aircraft ode by 0.25 sec (the result is stored in x)  */
double mg = 60000.0*9.81;
double mi = 1.0/60000;
auto aircraft_post = [] (state_type &x, const input_type &u) {
        /* the ode describing the aircraft */
        auto rhs =[] (state_type& xx,  const state_type &x, const input_type &u) {
                double c=(1.25+4.2*u[1]);
                xx[0] = mi*(u[0]*std::cos(u[1])-(2.7+3.08*c*c)*x[0]*x[0]-mg*std::sin(x[1]));
                xx[1] = (1.0/(60000*x[0]))*(u[0]*std::sin(u[1])+68.6*c*x[0]*x[0]-mg*std::cos(x[1]));
                xx[2] = x[0]*std::sin(x[1]);
        };

  ode_solver(rhs,x,u,sDIM,tau,10);
};

/* we integrate the growth bound by 0.25 sec (the result is stored in r)  */
/* lischitz matrix */
double L[3][2];
//state_type v={{2.77*0.01,5.16*0.0001,7.4*0.001}};
state_type w={{.108,0.002,0}};
auto radius_post = [] (state_type &r, const state_type&, const input_type &u) {

        L[0][0]=-0.0019*(2.7+3.08*(1.25+4.2*u[1])*(1.25+4.2*u[1]));
        L[0][1]=9.81;

        L[1][0]=0.00447+0.00481*u[1];
        L[1][1]=0.005163;

        L[2][0]=.03;
        L[2][1]=83;

        /* the ode for the growth bound */
        auto rhs =[] (state_type& rr,  const state_type &r, const input_type &) { 
                rr[0] = L[0][0]*r[0]+L[0][1]*r[1]+w[0]; /* L[0][2]=0 */
                rr[1] = L[1][0]*r[0]+L[1][1]*r[1]+w[1]; /* L[1][2]=0 */
                rr[2] = L[2][0]*r[0]+L[2][1]*r[1]+w[2]; /* L[2][2]=0 */
        };

        ode_solver(rhs,r,u,sDIM,tau,10);

};

int main() {
        /* to measure time */
        TicToc tt;

        /****************************************************************************/
        /* construct grid for the state space */
        /****************************************************************************/
        /* setup the workspace of the synthesis problem and the uniform grid */
        /* grid node distance diameter */
        //state_type eta={{25.0/100,3*M_PI/180/100,56.0/100}};
        double c=1;
        state_type eta={{c*25.0/362,c*3*M_PI/180/66,c*56.0/334}};
        /* lower bounds of the hyper rectangle */
        state_type lb={{58,-3*M_PI/180,0}};
        /* upper bounds of the hyper rectangle */
        state_type ub={{83,0,56}}; // changes the 56 to 28 for debugging reason and making the example smaller
        /* measurement disturbances  */
        state_type z={{0.0125,0.0025/180*M_PI,0.05}};
        //for(size_t i=0; i< sDIM; i++) {
        //  lb[i]=lb[i]+eta[i]/2+z[i];
        //  ub[i]=ub[i]-eta[i]/2-z[i];
        //}
        scots::UniformGrid<state_type> ss(sDIM,lb,ub,eta,z);
        //scots::UniformGrid<state_type> ss(sDIM,lb,ub,eta);
        std::cout << "Unfiorm grid details:" << std::endl;
        std::cout << "eta factor C is:" << c << std::endl;
        ss.printInfo(1);

        /****************************************************************************/
        /* construct grid for the input space */
        /****************************************************************************/
        /* lower bounds of the hyper rectangle */
        input_type ilb={{0,0}};
        /* upper bounds of the hyper rectangle */
        input_type iub={{32000,8*M_PI/180}};
        /* grid node distance diameter */
        input_type ieta={{32000,8.0/9.0*M_PI/180}};
        scots::UniformGrid<input_type> is(iDIM,ilb,iub,ieta);
        is.printInfo(1);

        /* transition system to be computed */
        scots::TransitionSystem ts;


        std::cout << "Number of transitions: " << ts.getNoTransitions() << std::endl;

        tt.tic();
        std::cout << "scots::AbstractionGB started & Compute transition relation" << std::endl;
        scots::AbstractionGB<state_type,input_type> abs(ss,is,ts);
        tt.toc();
        std::cout << "scots::AbstractionGB ended \n"<< std::endl;





        tt.tic();
        std::cout << "abs.computeTransitionRelation started \n \n  " << std::endl;
        abs.computeTransitionRelation(aircraft_post, radius_post);
        std::cout << "Number of transitions: " << ts.getNoTransitions() << std::endl;
        tt.toc();
        std::cout << "\n abs.computeTransitionRelation ended \n" << std::endl;


//  /* define function to check if the cell is in the  target set?  */
        state_type x;
        auto target = [&](const size_t idx)->bool {
                ss.itox(idx,x);
                /* function returns 1 if cell associated with x is in target set  */
                if (  63 <= (x[0]-eta[0]/2.0) && (x[0]+eta[0]/2.0)<= 75
                      && -3*M_PI/180 <= (x[1]-eta[1]/2.0) &&  (x[1]+eta[1]/2.0) <= 0
                      && 0 <= (x[2]-eta[2]/2.0) &&  (x[2]+eta[2]/2.0) <= 2.5
                      &&  -0.91 <=  ( (x[0]+eta[0]/2.0) * std::sin(x[1]-eta[1]/2.0) )
                      )
                        return true;
                  return false;
        };

//        ss.clearAbstractSet();
//        ss.addIndices(target);
//        scots::IO::writeToFile(&ss,"target.scs");
//
//        ss.fillAbstractSet();
//        ss.remIndices(target);
//        ss.remGridPoints(overflow);
//        scots::IO::writeToFile(&ss,"problemdomain.scs");

        /////////////////////////////  S   O   L    V   E  GAME  //////////////////////////////////////////////////

        tt.tic();

        std::cout << "Solve game " << std::endl;

        scots::ReachabilityGame reach(ts);
        reach.solve(target);
        tt.toc();

        std::cout << "Size: " << reach.size() << std::endl;


        return 1;
}

