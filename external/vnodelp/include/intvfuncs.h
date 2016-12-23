/*151:*/
#line 607 "./linalg.w"

#ifndef INTVFUNC_H
#define INTVFUNC_H

#include "vector_matrix.h"


namespace v_bias{
/*142:*/
#line 421 "./linalg.w"

inline double rad(const interval&a)
{
return 0.5*v_bias::width(a);
}

/*:142*/
#line 615 "./linalg.w"

}
namespace v_blas{

using namespace v_bias;

/*140:*/
#line 400 "./linalg.w"

inline bool interior(const iVector&a,const iVector&b)
{
return equal(a.begin(),a.end(),b.begin(),v_bias::interior);
}

/*:140*/
#line 621 "./linalg.w"

/*139:*/
#line 386 "./linalg.w"

inline bool subseteq(const iVector&a,const iVector&b)
{
return equal(a.begin(),a.end(),b.begin(),v_bias::subseteq);
}




/*:139*//*141:*/
#line 411 "./linalg.w"

inline bool disjoint(const iVector&a,const iVector&b)
{
return equal(a.begin(),a.end(),b.begin(),v_bias::disjoint);
}


/*:141*/
#line 622 "./linalg.w"

/*144:*/
#line 438 "./linalg.w"

inline void width(pVector&r,const iVector&a)
{
transform(a.begin(),a.end(),r.begin(),v_bias::width);
}



/*:144*/
#line 623 "./linalg.w"

/*145:*/
#line 447 "./linalg.w"

inline void midpoint(pVector&r,const iVector&a)
{
transform(a.begin(),a.end(),r.begin(),v_bias::midpoint);
}

/*:145*/
#line 624 "./linalg.w"

/*146:*/
#line 454 "./linalg.w"

inline void midpoint(pMatrix&R,const iMatrix&A)
{
for(unsigned int i= 0;i<A.size();i++)midpoint(R[i],A[i]);
}





/*:146*/
#line 625 "./linalg.w"

/*147:*/
#line 471 "./linalg.w"

inline bool intersect(iVector&z,const iVector&x,const iVector&y)
{
interval c;
for(unsigned int i= 0;i<sizeV(y);i++)
{
bool b= v_bias::intersect(c,x[i],y[i]);
if(!b)return false;
else z[i]= c;
}
return true;
}




/*:147*/
#line 626 "./linalg.w"

/*143:*/
#line 429 "./linalg.w"

inline void rad(pVector&r,const iVector&v)
{
transform(v.begin(),v.end(),r.begin(),v_bias::rad);
}


/*:143*/
#line 627 "./linalg.w"

double compH(const iVector&a,
const iVector&b);
}

#endif


/*:151*/
