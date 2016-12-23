
#include <ostream> 
#include "vnode.h"

using namespace std;
using namespace vnodelp;

template<typename var_type> 
void R3body(int n, var_type*Yp, const var_type*Y, var_type t,
	    void*param)
{
  // We encode as close as possible to the AWA encoding 
#define U1   Y[0]
#define U2   Y[1]
#define U3   Y[2]
#define U4   Y[3]
  
#define F1   Yp[0]
#define F2   Yp[1]
#define F3   Yp[2]
#define F4   Yp[3]
  
  var_type x = U1,  y = U2,  x1 = U3,  y1=U4;
  
  interval m  = 1.0/string_to_interval("82.45");
  interval l  = 1.0 - m;
  
  var_type xm = x+m;
  var_type xl = x-l;
  
  var_type ys = sqr(y);
  var_type rm = pow( sqr(xm) + ys, interval(-1.5) );
  var_type rl = pow( sqr(xl) + ys, interval(-1.5) );

  F1 = x1;
  F2 = y1;
  F3 = x + 2.0*y1 - l*xm*rm - m*xl*rl;
  F4 = y - 2.0*x1 - l*y *rm - m*y *rl;
}


int main()
{
  const int n = 4;
  iVector y(n);
  
  y[0] = string_to_interval("1.2");
  y[1] = 0;
  y[2] = 0;
  y[3] = string_to_interval("-1.04935750983");
  
  AD *ad = new FADBAD_AD(n,R3body,R3body);
  VNODE *Solver = new VNODE(ad);
  
  Solver->setTols(1e-15,1e-15);
  
  interval t = 0.0, tend = string_to_interval("6.192169331396");
  Solver->integrate(t,y,tend);
  
  // enclosure computed by AWA
  iVector awa(n);
  double l = inf(string_to_interval("1.19999999999437E+000"));
  double r = sup(string_to_interval("1.20000000000552E+000"));
  awa[0] = interval(l,r);
  
  l = inf(string_to_interval("-9.30335015353590E-011"));
  r = sup(string_to_interval("-6.79911866837946E-011"));
  awa[1] = interval(l,r);
  
  l = inf(string_to_interval("-1.51974042204144E-010"));
  r = sup(string_to_interval("-1.28985217342847E-010"));
  awa[2] = interval(l,r);
  
  l = inf(string_to_interval("-1.04935750984473E+000"));
  r = sup(string_to_interval("-1.04935750981524E+000"));
  awa[3] = interval(l,r);
  
  
  // y and awa must intersect
  if ( !intersect(y, awa, y) )
    {
      cerr << "Test n-5 failed: y and y1 must intersect" << endl;
      cout << "t = " << t << endl;
      printVector(awa, "awa ");
      printVector(y,  "y  ");
      return -1;
    }
  
  cout << "Test n-5 SUCCESSFUL" << endl;
  return 0;
}




