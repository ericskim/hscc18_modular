/*349:*/
#line 148 "./fadbad.w"

#ifndef Fadbad_ODE
#define Fadbad_ODE



#include "vnodeinterval.h"
#include "basiclinalg.h"
#include "ad_ode.h"

#include "ffadiff.h"
#include "fadbad_intv.inc"



namespace vnodelp{
/*339:*/
#line 11 "./fadbad.w"



typedef T<interval> Tinterval;
typedef void(*Tfunction)(int n,
Tinterval*yp,
const Tinterval*y,
Tinterval t,
void*param);


class FadbadODE:public AD_ODE
{
public:

FadbadODE(int n,Tfunction,void*param= 0);

void set(const interval&t0,const iVector&y0,
const interval&h,int k);
void compTerms();
void sumTerms(iVector&sum,int m);
void getTerm(iVector&term,int i)const;
interval getStepsize()const;

void eval(void*param);

~FadbadODE();
private:
Tfunction fcn;
Tinterval*y_coeff,*f_coeff,t;

int size;
int order;
interval stepsize;
};


/*:339*/
#line 164 "./fadbad.w"

}
#endif

/*:349*/
