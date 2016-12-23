/*195:*/
#line 181 "./adabst.w"

#ifndef AD_VAR_H
#define AD_VAR_H
#include "vnodeinterval.h"
#include "vector_matrix.h"
namespace vnodelp{
using namespace v_bias;
using namespace v_blas;
/*189:*/
#line 77 "./adabst.w"

class AD_VAR
{
public:
virtual void set(const v_bias::interval&t0,
const iVector&y0,const v_bias::interval&h,int k)= 0;

virtual void compTerms()= 0;
virtual void sumTerms(iMatrix&sum,int m)= 0;
virtual void getTerm(iMatrix&term,int i)const= 0;

virtual void eval(void*param)= 0;
virtual~AD_VAR(){}
};


/*:189*/
#line 189 "./adabst.w"

}
#endif

/*:195*/
