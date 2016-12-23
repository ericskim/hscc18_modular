/*374:*/
#line 120 "./misc.w"

#ifndef MISCFUNS_H
#define MISCFUNS_H
#include <cmath> 

#include "basiclinalg.h"
using namespace std;

/*367:*/
#line 21 "./misc.w"

namespace v_bias{
inline bool finite_interval(const interval&a)
{
return(isfinite(inf(a))&&isfinite(sup(a)));
}
}

/*:367*//*368:*/
#line 30 "./misc.w"

namespace v_blas{
inline bool finite_interval(const iVector&a)
{
for(unsigned int i= 0;i<a.size();i++)
if(!v_bias::finite_interval(a[i]))return false;
return true;
}
}

/*:368*/
#line 128 "./misc.w"



namespace vnodelp{
void vnodeMessage(const char*s);
}

#endif





#line 1 "./timing.w"


/*:374*/
