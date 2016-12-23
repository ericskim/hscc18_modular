/*194:*/
#line 167 "./adabst.w"


#ifndef AD_ODE_H
#define AD_ODE_H
#include "vnodeinterval.h"
#include "vector_matrix.h"
namespace vnodelp{
using namespace v_bias;
using namespace v_blas;
/*187:*/
#line 16 "./adabst.w"

class AD_ODE
{
public:
virtual void set(const v_bias::interval&t0,const iVector&y0,
const v_bias::interval&h,int k)= 0;

virtual void compTerms()= 0;
virtual void sumTerms(iVector&sum,int m)= 0;
virtual void getTerm(iVector&term,int i)const= 0;
virtual v_bias::interval getStepsize()const= 0;

virtual void eval(void*param)= 0;

virtual~AD_ODE(){}
};

/*:187*/
#line 176 "./adabst.w"

}
#endif

/*:194*/
