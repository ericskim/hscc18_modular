/*
 * compute_enclosure.cc
 *
 *  created: Dec 2016
 *   author: Matthias Rungger
 */

/* 
 * see readme 
 */

#include <ostream> 
#include "vnode.h"

const int state_dim = 3; 
const int input_dim = 2; 

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
  /* set measurement errors */
  vnodelp::iVector z(state_dim);
  z[0] = vnodelp::interval(-0.0125,.0125);
  z[1] = vnodelp::interval(-0.0025,0.0025)*pi()/180.0;
  z[2] = vnodelp::interval(-0.05,0.05);
  /* set state space bounds */
  vnodelp::iVector x(state_dim);
  vnodelp::interval p = -3.0*pi()/180.0;
  x[0] = vnodelp::interval(58.0,83.0)+z[0];
  x[1] = vnodelp::interval(inf(p),0)+z[1];
  x[2] = vnodelp::interval(0,56)+z[2];
  /* set input space bounds */
  param_t* param = new param_t;
  param->u1 = vnodelp::interval(0,32000);
  param->u2 = vnodelp::interval(0,sup(8.0*pi()/180.0));

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
  /* compute apriori enclosure */
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


  std::cout << "Apriori enclosure computed for " << std::endl;
  std::cout << "sampling interval " << interval(inf(t0),tau) << std::endl;
  std::cout << "input bounds " << param->u1 << "x" << param->u2 << std::endl;
  std::cout << "state bounds " << x[0] << "x" << x[1] << "x" << x[2] << std::endl;
  std::cout << std::endl;
  
  vnodelp::iVector y=hoe->getApriori();
  std::cout << "enclosure " << y[0] << "x" << y[1] << "x" << y[2] << std::endl;

  delete ad;
  delete control;
  delete hoe;
  delete param;

  return 0;
}
