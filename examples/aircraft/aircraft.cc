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

/* we integrate the aircraft ode by 0.25 sec (the result is stored in x)  */
double mg = 60000.0*9.81;
double mi = 1.0/60000;
auto aircraft_post = [] (state_type &x, input_type &u)->void {
        /* the ode describing the aircraft */
        auto rhs =[] (state_type& xx,  const state_type &x, input_type &u)->void {
                double c=(1.25+4.2*u[1]);
                xx[0] = mi*(u[0]*std::cos(u[1])-(2.7+3.08*c*c)*x[0]*x[0]-mg*std::sin(x[1]));
                xx[1] = (1.0/(60000*x[0]))*(u[0]*std::sin(u[1])+68.6*c*x[0]*x[0]-mg*std::cos(x[1]));
                xx[2] = x[0]*std::sin(x[1]);
        };

        size_t nint=5; /* number of intermediate step size 5*/
        double h=0.05; /* h* nint = sampling time 0.05*5=0.25*/
        ode_solver(rhs,x,u,nint,h);
};

/* we integrate the growth bound by 0.25 sec (the result is stored in r)  */
/* lischitz matrix */
double L[3][2];
//state_type v={{2.77*0.01,5.16*0.0001,7.4*0.001}};
state_type w={{.108,0.002,0}};
auto radius_post = [] (state_type &r, input_type &u)->void {

        L[0][0]=-0.0019*(2.7+3.08*(1.25+4.2*u[1])*(1.25+4.2*u[1]));
        L[0][1]=9.81;

        L[1][0]=0.00447+0.00481*u[1];
        L[1][1]=0.005163;

        L[2][0]=.03;
        L[2][1]=83;

        /* the ode for the growth bound */
        auto rhs =[] (state_type& rr,  const state_type &r, input_type u)->void { //w u l input type ??? added
                rr[0] = L[0][0]*r[0]+L[0][1]*r[1]+w[0]; /* L[0][2]=0 */
                rr[1] = L[1][0]*r[0]+L[1][1]*r[1]+w[1]; /* L[1][2]=0 */
                rr[2] = L[2][0]*r[0]+L[2][1]*r[1]+w[2]; /* L[2][2]=0 */
        };

        size_t nint=5; /* number of intermediate step size */
        double h=0.05; /* h* nint = sampling time */
        ode_solver(rhs,r,u,nint,h);

        //for(size_t i=0; i<sDIM; i++)
        //  r[i]+=v[i];

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
        scots::AbstractionGB<state_type,input_type> abs(&ss,&is,&ts);
        tt.toc();
        std::cout << "scots::AbstractionGB ended \n"<< std::endl;


        auto overflow = [] (state_type)->bool { return false; };



        tt.tic();
        std::cout << "abs.computeTransitionRelation started \n \n  " << std::endl;
        abs.computeTransitionRelation(aircraft_post, radius_post, overflow,"preOnly" ); //"postOnly" "preOnly"
        std::cout << "Number of transitions: " << ts.getNoTransitions() << std::endl;
        tt.toc();
        std::cout << "\n abs.computeTransitionRelation ended \n" << std::endl;


// scots::IO::writeToFile(&ts,"ts.scs"); // take one hour !


//begin by Meysam commented :
//  size_t no=0;
// for(size_t i=0; i < ts.getN() ;i++)
//  for(size_t j=0; j < ts.getM() ;j++)
//   no = ( no < ts.relation_[i].label[j].no ? ts.relation_[i].label[j].no : no );

//  std::cout << "Avg no of post: " << no << std::endl;

// end by Meysam commented




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

        /////////////////////////////  S   O   L    V   E  GAME  //////////////////////////////////////////////////

        tt.tic();

        std::cout << "Solve game " << std::endl;

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
                for(size_t i=0; i<sDIM; i++)
                        tmp[i]=x[i]+h/2*k[0][i];

                rhs(k[1],tmp, u);
                for(size_t i=0; i<sDIM; i++)
                        tmp[i]=x[i]+h/2*k[1][i];

                rhs(k[2],tmp, u);
                for(size_t i=0; i<sDIM; i++)
                        tmp[i]=x[i]+h*k[2][i];

                rhs(k[3],tmp, u);
                for(size_t i=0; i<sDIM; i++)
                        x[i] = x[i] + (h/6)*(k[0][i] + 2*k[1][i] + 2*k[2][i] + k[3][i]);
        }
}
