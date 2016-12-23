/*360:*/
#line 334 "./fadbad.w"

#ifndef Fadbad_Var_ODE
#define Fadbad_Var_ODE



#include "vnodeinterval.h"
#include "vector_matrix.h"
#include "ad_var.h"



#include "ffadiff.h"
#include "fadbad_intv.inc"



namespace vnodelp{
using namespace v_bias;
using namespace v_blas;
/*352:*/
#line 194 "./fadbad.w"


typedef T<F<interval> > TFinterval;

typedef void(*TFfunction)(int n,TFinterval*yp,const TFinterval*y,
TFinterval tf,
void*param);



class FadbadVarODE:public AD_VAR
{
public:

FadbadVarODE(int n,TFfunction f,void*param= 0);

void set(const interval&t0,const
iVector&y0,const interval&h,int k);
void compTerms();
void sumTerms(iMatrix&sum,int m);
void getTerm(iMatrix&term,int i)const;
void eval(void*param){fcn(size,tf_out,tf_in,tf,param);}
~FadbadVarODE();
private:

TFinterval*tf_in,*tf_out,tf;

TFfunction fcn;
int size;
int order;
interval stepsize;
};

/*:352*/
#line 354 "./fadbad.w"
}
#endif

/*:360*/
