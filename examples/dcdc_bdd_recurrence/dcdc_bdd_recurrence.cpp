#include <iostream>
#include <array>
#include <cmath>

/* CUDD header */
#include "cudd.h"
#include "cuddObj.hh"

/* SCOTS_BDD define */
#ifndef SCOTS_BDD
#define SCOTS_BDD
#endif

/* SCOTS header */
#include "scots.hh"

/* ode solver */
#include "RungeKutta4.hh"

/* time profiling */
#include "TicToc.hh"

/* memory profiling */
#include <sys/time.h>
#include <sys/resource.h>
struct rusage usage;


/***********************************************
 *
 * define state space and dynamic
 *
 * ********************************************/

/* state space dim */
const int state_dim=2;
/* input space dim */
const int input_dim=1;
/* sampling time */
const double tau = 0.5;

/*
 * data types for the elements of the state space
 * and input space used by the ODE solver
 */
using state_type = std::array<double,state_dim>;
using input_type = std::array<double,input_dim>;

/* parameters for system dynamics */
const double xc=70;
const double xl=3;
const double rc=0.005;
const double rl=0.05;
const double ro=1;
const double vs=1;

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
    scots::runge_kutta_fixed4(rhs,x,u,state_dim,tau,5);
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
    scots::runge_kutta_fixed4(rhs,r,u,state_dim,tau,5);
};

int main()
{

    std::cout << "dcdc example for the use of scots with bdd: \n"
              << "to solve a recurrence problem for the dcdc converter, \n"
              << "we we implement the nested fixed point algorithm: \n"
              << "nu Z.intersection_{i=1-N}(mu Y[i].((T[i] & pre(Z)) | pre(Y[i])) \n";

    std::cout << "\n press Enter to show uniform grids and Targets: \n";
    std::getchar();

    /* to measure time */
    TicToc tt;
    /* BDD manager */
    Cudd manager;
    /* enable variable reordering */
    manager.AutodynEnable();

     /***********************************************
     *
     * setup of the uniform grids
     *
     * ********************************************/

    /* construct grid for the state alphabet */
    /* grid node distance diameter */
    state_type eta={{20.0/4e3,20.0/4e3}};
    /* lower bounds of the hyper-rectangle */
    state_type lb={{0.649,4.949}};
    /* upper bounds of the hyper-rectangle */
    state_type ub={{1.65,5.95}};
    scots::SymbolicSet ss_pre(manager,state_dim,lb,ub,eta);
    scots::SymbolicSet ss_post(manager,state_dim,lb,ub,eta);
    std::cout << "\n state alphabet details: \n";
    ss_pre.print_info();

    /* construct grid for the input alphabet */
    /* hyper-rectangle [1,2] with grid node distance 1 */
    scots::SymbolicSet ss_input(manager,input_dim,input_type{{.99}},input_type{{2.1}},input_type{{1}});
    std::cout << "\n input alphabet details: \n";
    ss_input.print_info();

    /***********************************************
     *
     * Definition of the Targets:
     * Targets are represented as rectangles in the
     * state space and defined by lambda functions
     *
     * we create N=2 targets
     *
     * ********************************************/

    /* Vectors for the Targets and their individual controller */
    std::vector<BDD> targets;
    std::vector<BDD> controller;

    /* lambda for target 1*/
    auto target_1 = [&ss_pre,&eta](const scots::abs_type& idx) {
        /* bounding rectangle for T_1*/
        double h[4] = {1.3,1.6,5.7,5.9};
        state_type x;
        ss_pre.itox(idx,x);
        double c1= eta[0]/2.0+1e-10;
        double c2= eta[1]/2.0+1e-10;
        if ((h[0]+c1) <= x[0] && x[0] <= (h[1]-c1) &&
                (h[2]+c2) <= x[1] && x[1] <= (h[3]-c2)) {
            return true;
        }
        return false;
    };
    /* lambda for target 2*/
    auto target_2 = [&ss_pre,&eta](const scots::abs_type& idx) {
        /* bounding rectangle for T_2*/
        double h[4] = {0.8,1.1,5.0,5.2};
        state_type x;
        ss_pre.itox(idx,x);
        double c1= eta[0]/2.0+1e-10;
        double c2= eta[1]/2.0+1e-10;
        if ((h[0]+c1) <= x[0] && x[0] <= (h[1]-c1) &&
                (h[2]+c2) <= x[1] && x[1] <= (h[3]-c2)) {
            return true;
        }
        return false;
    };
    /* fill the target set with T_1 - T_N */
    targets.push_back(ss_pre.ap_to_bdd(manager,target_1));
    targets.push_back(ss_pre.ap_to_bdd(manager,target_2));

    /* initiate N=2 controller as empty */
    controller.push_back(manager.bddZero());
    controller.push_back(manager.bddZero());

    std::cout << "\n targets: \n"
              << "Target1: [1.3,1.6];[5.7,5.9] \n"
              << "Target2: [0.8,1.1];[5.0,5.2] \n";


    /***********************************************
     *
     * Compute transition function of symbolic model
     *
     * ********************************************/
    std::cout << "\n press Enter to compute the transition function: \n";
    std::getchar();
    std::cout << "computing the transition function:\n";

    tt.tic();
    /* SymbolicModel class to compute the BDD encoding the transition function */
    scots::SymbolicModel<state_type,input_type> sym_model(ss_pre,ss_input,ss_post);
    /* flag for the number of transitions*/
    size_t no_trans;
    /* computing the growth bound*/
    BDD TF = sym_model.compute_gb(manager,system_post,radius_post,no_trans);
    tt.toc();

    std::cout << "number of transitions " << no_trans  << "\n";
    if(!getrusage(RUSAGE_SELF, &usage)) {
        std::cout << "memory per transition: " << usage.ru_maxrss/(double)no_trans<< "\n";
    }

    manager.DebugCheck();

    /***********************************************
     *
     * Synthesis
     *
     * ********************************************/
    std::cout << "\n press Enter to continue with the synthesis: ";
    std::getchar();

    /* class for extracting the pre lists*/
    scots::EnfPre enf_pre(manager,TF,sym_model);
    size_t i;

    /* outer fp*/
    BDD Z   = manager.bddZero();
    BDD ZZ  = manager.bddOne();

    /* inner fps*/
    std::vector<BDD> Y;
    std::vector<BDD> YY;
    for(i = 0; i < 2; i++)
    {
        Y.push_back(manager.bddOne());
        YY.push_back(manager.bddZero());
    }

    /* helper */
    BDD U=ss_input.get_cube(manager);

    /***********************************************
     *
     * Fixed Point Computation:
     *
     * nu Z.intersection_{i=1-N}(mu Y[i].((T[i] & pre(Z)) | pre(Y[i]))
     *
     * ********************************************/
    tt.tic();
    while(Z != ZZ)
    {
        Z = ZZ;
        BDD preZ = enf_pre(Z);
        ZZ = manager.bddOne();
        for(int i = 0; i < 2; i++)
        {
            YY[i] = manager.bddZero();
            controller[i] = manager.bddZero();
            while(Y[i] != YY[i])
            {
                Y[i]   = YY[i];
                YY[i]  = (targets[i] & preZ) | (enf_pre(Y[i]));

                BDD N = YY[i] & (!(controller[i].ExistAbstract(U)));
                controller[i] = controller[i] | N;
            }

            ZZ = ZZ & YY[i];
        }
    }
    tt.toc();

    /*checking on the sizes of the winning domains*/
    for(i = 0; i < 2; i++)
    {
        std::cout << "winning domain size of Target" + std::to_string(i+1) + " : "
                  << ss_pre.get_size(manager,controller[i]) << "\n";
    }

    /***********************************************
     *
     * Saving of the controller
     *
     * ********************************************/
    std::cout << "\n press Enter to save the controllers: ";
    std::getchar();

    scots::SymbolicSet symbSet(ss_pre,ss_input);
    for(i = 0; i < 2; i++)
    {
        std::cout << "write controller" + std::to_string(i+1) + " to file : ";
        if(write_to_file(manager,symbSet,controller[i],"controller"+std::to_string(i+1)))
        {
            std::cout << "done. \n";
        }
        else
        {
            std::cout << "failed. \n";
        }
    }

    /***********************************************
     *
     * Simulation
     *
     * ********************************************/
    std::cout << "\n simulation:";

    /* input vectors*/
    input_type u ={{0}};
    std::vector<double> inputs;

    /* Initial state */
    state_type x={{0.7, 5.4}};

    std::cout << "\n press Enter to start the simulation to Target1: [1.3,1.6];[5.7,5.9] \n";
    std::getchar();

    while(1)  {
        std::cout << x[0] <<  " "  << x[1] << "\n";
        /*calculate the abstract index and check if we are in the target set of Target1*/
        scots::abs_type idx = ss_pre.xtoi(x);
        if(target_1(idx))
        {
            std::cout << " \n -> reached Target1: [1.3,1.6];[5.7,5.9]  \n \n";
            break;
        }
        /* returns a std vector with the valid control inputs */
        inputs = (symbSet.restriction(manager,controller[0],x));

        if(!inputs.empty())
        {
            /*we always take the first in his example*/
            u[0] = inputs[0];
            /*simulate to the next state*/
            system_post(x,u);
        }
        else
        {   /* we are out of the winning domain*/
            std::cout << " \n ->error: out of Domain \n";
            return 0;
        }
    }

    std::cout << "\n press Enter to continue the simulation to Target2: [0.8,1.1];[5.0,5.2] \n";
    std::getchar();

    while(1) {
        std::cout << x[0] <<  " "  << x[1] << "\n";
        /*calculate the abstract index and check if we are in the target set of Target2*/
        scots::abs_type idx = ss_pre.xtoi(x);
        if(target_2(idx))
        {
            std::cout << "\n -> reached Target2: [0.8,1.1];[5.0,5.2]  \n \n";
            break;
        }
        /* returns a std vector with the valid control inputs */
        inputs = (symbSet.restriction(manager,controller[1],x));

        if(!inputs.empty())
        {
            /*we always take the first in his example*/
            u[0] = inputs[0];
            /*simulate to the next state*/
            system_post(x,u);
        }
        else
        {   /* we are out of the winning domain*/
            std::cout << " \n ->error: out of Domain \n";
            return 0;
        }
    }

    return 1;
}


