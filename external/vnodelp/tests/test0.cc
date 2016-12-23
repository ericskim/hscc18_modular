#include <ostream> 
#include "vnode.h"

using namespace std;
using namespace vnodelp;

/* Make sure the functions VNODE-LP uses compile and execute. */
int main()
{
  interval x(10.0,10.1);
  interval y;

  double d;
  bool b;
  
  b = inf(x);
  b = sup(x);
  
  d = midpoint(x);
  d = width(x);
  d = mag(x);
  
  b = subseteq (x,y);
  b = interior (x,y);
  b = disjoint (x,y);
  b = intersect(x,x,x);
  
  y = pi();
  
  y = pow(x,y);
  y = pow(x,10);
  
  y = exp(x);
  y = log(x);
  
  y = sqr(x);
  y = sqrt(x);
  
  y = sin(x);
  y = cos(x);
  y = tan(x);
  
  x = interval(-0.1,0.1);
  y = asin(x);
  y = acos(x);
  y = atan(x);
  
  cout << "Test 0 SUCCESSFUL" << endl;
  return 0;
}
