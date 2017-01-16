/*
 * compute_enclosure_per_input.cc
 *
 *  created on: 07.12.2016
 *  author: rungger
 */

/* 
 * see readme 
 */
#include <iostream> 
#include <array>

#include "vnode.h"
#include "scots.hh"

const int state_dim = 3; 
const int input_dim = 2; 

typedef std::array<double,input_dim> input_type;

struct param_t {
  vnodelp::interval u1;
  vnodelp::interval u2;
};

template<typename var_type> 
void aircraft(int n, var_type* xx, const var_type* x, var_type t, void* param) {

  param_t* p = (param_t*) param;
  vnodelp::iVector u(input_dim);
  u[0]= p->u1;
  u[1]= p->u2;

  vnodelp::interval mg = vnodelp::interval(60000.0)*9.81;
  vnodelp::interval mi = vnodelp::interval(1.0)/60000.0;
  vnodelp::interval c = vnodelp::interval(1.25+4.2*u[1]);

  vnodelp::interval w0 = vnodelp::interval(-0.108,0.108);
  vnodelp::interval w1 = vnodelp::interval(-0.002,0.002);

  xx[0] = mi*(u[0]*cos(u[1])-(2.7+3.08*c*c)*x[0]*x[0]-mg*sin(x[1]))+w0;
  xx[1] = (1.0/(60000.0*x[0]))*(u[0]*sin(u[1])+68.6*c*x[0]*x[0]-mg*cos(x[1]))+w1;
  xx[2] = x[0]*sin(x[1]);
}

int main(){

  /* sampling time */
  double tau = 0.25;
  /* set state space bounds */
  vnodelp::iVector x(state_dim);
  vnodelp::interval p = -3.0*pi()/180.0;
  x[0] = vnodelp::interval(58.0,83.0);
  x[1] = vnodelp::interval(inf(p),0);
  x[2] = vnodelp::interval(0,56);
  /* set input space bounds */
  /* lower bounds of the hyper rectangle */
  input_type ilb={{0,0}};
  /* upper bounds of the hyper rectangle */
  input_type iub={{32000,8*M_PI/180}};
  /* grid node distance diameter */
  input_type ieta={{32000,8.0/9.0*M_PI/180}};
  scots::UniformGrid is(input_dim,ilb,iub,ieta);
  /* setup param struct to update input */
  param_t* param = new param_t;

  /* set automated differentiation object */
  vnodelp::AD* ad= new vnodelp::FADBAD_AD(state_dim,aircraft,aircraft,param);
  /* setup control obj */
  vnodelp::Control* control = new vnodelp::Control;
  /* setup computation obj for apriori enclosure */
  vnodelp::HOE* hoe = new vnodelp::HOE(state_dim);
  hoe->setTrialStepsize(tau);
  hoe->setTrialOrder(control->order);
  hoe->set(control,ad);

  /* init apriori computation */
  vnodelp::interval t0 = 0;
  hoe->init(t0,x);

  std::cout << "Apriori enclosure computation for " << std::endl;
  std::cout << "sampling interval " << interval(inf(t0),tau) << std::endl;
  std::cout << "state bounds " << x[0] << "x" << x[1] << "x" << x[2] << std::endl;
  std::cout << std::endl;

  input_type u;

  /* loop over all inputs */
  for(size_t i=0; i<is.size(); i++) {

    /****************************************/
    /* first compute a priori enclosure     */
    /****************************************/

    /* get input */
    is.itox(i,u);
    param->u1 = vnodelp::interval(u[0]);
    param->u2 = vnodelp::interval(u[1]);
    std::cout << "input: " << param->u1 << "x" << param->u2 << std::endl;

    /* tell solver that params have changed */
    ad->eval(param); 

    /* compute apriori enclosure */
    hoe->setTrialStepsize(tau);
    bool info;
    hoe->compAprioriEnclosure(t0,x,info);
    if(info==false) {
      std::cout << "Could not validate solution" << std::endl;
      return 0;
    }
    /* check if apriori enclosure was successfully computed for [0,tau] */
    vnodelp::interval Tj= hoe->getTrialT();
    if(v_bias::subseteq(vnodelp::interval(inf(t0),tau),Tj)) {
      hoe->acceptSolution();
    } else {
      std::cout << "Could not validate solution" << std::endl;
      return 0;
    }

    vnodelp::iVector y=hoe->getApriori();
    std::cout << "enclosure " << std::endl;
    std::cout << y[0] << std::endl;
    std::cout << y[1] << std::endl;
    std::cout << y[2] << std::endl;


    /****************************************/
    /* compute growth bound by bounding */
    /****************************************/
    /* setup computation of jacobian */
    vnodelp::iMatrix M;
    vnodelp::sizeM(M,state_dim);

    ad->tayl_coeff_var->set(1.0,y,1.0,1); 
    ad->tayl_coeff_var->compTerms();
    ad->tayl_coeff_var->getTerm(M,1); 

    std::cout << "growth bound " << std::endl;
    for(int i=0; i<state_dim; i++) {
      for(int j=0; j<state_dim; j++) {

        if(i==j)
          std::cout << interval(sup(M[i][j])) <<  "       ";
        else
          std::cout <<interval(mag(M[i][j]))<< "       ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << std::endl;


  }
  delete ad;
  delete control;
  delete hoe;
  delete param;

  return 0;
}
