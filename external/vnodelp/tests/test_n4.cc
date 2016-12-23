#include <ostream> 
#include "vnode.h"

using namespace std;
using namespace vnodelp;

template<typename var_type> 
void Pend(int n, var_type*yp, const var_type*y, var_type t,
	  void*param)
{
  double b = 0;
  double c = 0;
  
  yp[0] = y[1];
  yp[1] = b*cos(t) - c*y[1] - sin(y[0]);
}


int main()
{
  const int n = 2;
  iVector y(n);
  
  y[0] = 0.0;
  double l = inf(string_to_interval("1.9999"));
  double r = sup(string_to_interval("2.0001"));
  y[1] = interval(l,r);
  
  AD *ad = new FADBAD_AD(n,Pend,Pend);
  VNODE *Solver = new VNODE(ad);
  
  Solver->setTols(1e-15,1e-15);
  
  interval t = 0.0, tend = 8.0;
  Solver->integrate(t,y,tend);
  
  iVector awa(n);
  
  // These are the intervals computed by AWA.
  l = inf(string_to_interval("3.06568304887608E+000"));
  r = sup(string_to_interval("3.21481855738090E+000"));
  awa[0] = interval(l,r);
  
  l = inf(string_to_interval("-7.32467052611495E-002"));
  r = sup(string_to_interval("7.59304059823085E-002"));
  awa[1] = interval(l,r);
  
  // y and awa must intersect
  if ( !intersect(y, awa, y) )
    {
      cerr << "Test n-4 failed: y and y1 must intersect"<< endl;
      printVector(awa, "awa ");
      printVector(y,  "y  ");
      return 1;
    }
  
  cout << "Test n-4 SUCCESSFUL" << endl;
  return 0;
}

