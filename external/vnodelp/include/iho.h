/*292:*/
#line 1625 "./iho.w"


#ifndef IHO_H
#define IHO_H
#include <cassert> 
#include  <cmath> 
#include "ad_ode.h"
#include "ad_var.h"
#include "hoe.h"
#include "intvfuncs.h"
#include "basiclinalg.h"
#include "matrixinverse.h"

namespace vnodelp{
using namespace v_bias;
using namespace v_blas;
/*235:*/
#line 273 "./iho.w"


class IHO
{
public:

IHO(int n);
/*277:*/
#line 1458 "./iho.w"


void set(Control*c,HOE*hoem,AD*ad0)
{
assert(c&&hoem&&ad0);
control= c;
hoe= hoem;
ad= ad0;
}


void getTightEnclosure(iVector&y)const{y= solution->y;}

void init(const interval&t,const iVector&y)
{
trial_solution->init(t,y);
solution->init(t,y);
}

const iVector&getGlobalExcess()const{return globalExcess;}



/*:277*/
#line 280 "./iho.w"


void compCoeffs();

void compTightEnclosure(interval&t_next);

void acceptSolution();


virtual~IHO();

private:
void compCpq(int p,int q);
void compCqp(int p,int q);
interval compErrorConstant(int p,int q);

private:
unsigned int p,q,order_trial;
interval h_trial;

pVector y_pred_point;

iVector y,y_pred,globalExcess,temp,temp2,x,u_next,predictor_excess,
corrector_excess,z,w,gj,term,d,s;

iMatrix Fj,M,Cinv,G,B,S,A,Q,U,V,Ainv;
pMatrix C,A_point;


MatrixInverse*matrix_inverse;
Solution*solution,*trial_solution;
interval*C_pq,*C_qp;

HOE*hoe;
AD*ad;
Control*control;
interval errorConstant;
};






/*:235*/
#line 1641 "./iho.w"


}

#endif



/*:292*/
