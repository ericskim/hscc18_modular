/*152:*/
#line 636 "./linalg.w"

#include <limits> 
#include <algorithm> 
#include "vnodeinterval.h"
#include "vnoderound.h"
#include "vector_matrix.h"




namespace v_blas{
/*150:*/
#line 584 "./linalg.w"

/*149:*/
#line 553 "./linalg.w"

#include <climits> 
using namespace std;
using namespace v_bias;
inline double compH(const v_bias::interval&a,const v_bias::interval&b)
{
if(inf(a)==0&&sup(a)==0)
return numeric_limits<double> ::max();

round_down();
if(inf(a)>=0)
return sup(b)/sup(a);

return inf(b)/inf(a);
}




/*:149*/
#line 585 "./linalg.w"


double compH(const iVector&a,
const iVector&b)
{
double hmin= compH(a[0],b[0]);
for(unsigned int i= 1;i<sizeV(a);i++)
{
double h= compH(a[i],b[i]);

if(h<hmin)
hmin= h;
}
return hmin;

}




/*:150*/
#line 647 "./linalg.w"

}
#line 1 "./qr.w"
/*:152*/
