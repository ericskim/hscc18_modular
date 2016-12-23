/*373:*/
#line 96 "./misc.w"

#include <ostream> 
#include <stdarg.h> 
#include <stdio.h> 

#include <stdlib.h> 
#include <iostream> 
using namespace std;

#include "vnodeinterval.h"
#include "basiclinalg.h"

using namespace std;
namespace vnodelp{
using namespace v_bias;
using namespace v_blas;
/*371:*/
#line 67 "./misc.w"


void checkIntersection(const iVector&a,const iVector&b)
{
v_bias::interval ai,bi;
for(unsigned int i= 0;i<vnodelp::sizeV(a);i++)
{
ai= a[i];
bi= b[i];
if(vnodelp::disjoint(ai,bi))
cout<<"i = "<<i<<"    "<<endl<<ai<<endl<<bi;
}
}

/*:371*/
#line 112 "./misc.w"

/*370:*/
#line 59 "./misc.w"

void vnodeMessage(const char*s)
{
cerr<<endl<<" *** VNODE-LP: "<<s<<endl;
}

/*:370*/
#line 113 "./misc.w"

}




/*:373*/
