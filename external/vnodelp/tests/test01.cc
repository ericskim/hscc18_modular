#include <ostream> 
#include "vnode.h"

using namespace std;
using namespace vnodelp;

/* Make sure the functions FADBAD++ uses compile and execute. */
template<typename var_type> 
void DummyODE(int n, var_type*yp, const var_type*y, var_type t,
	      void*param)
{
  yp[0] =  sqrt(y[0]);
  yp[1] = -sqr (y[1]);    
  
  yp[2] =  -exp(y[2]);
  yp[3] =  log (y[3]);
  
  yp[4] = -sin(y[4]);
  yp[5] = -asin(1.0/(t+10.0))*y[4];
  
  yp[5] = -cos(y[5]);
  yp[6] = -acos(1.0/(t+10.0))*y[6];
  
  yp[7] = -tan(y[7]);
  yp[8] = -atan(1.0/(t+10.0))*y[8];
  
  yp[9]  = -pow(y[9],1.0);
  yp[10] = -pow(y[10],interval(1.0,1.01));
}

int main()
{
  const int n = 11;
  
  interval t = 0.0, tend = 1;
  iVector y(n);
  setV(y,1.0);
  
  AD *ad = new FADBAD_AD(n,DummyODE,DummyODE);
  VNODE *Solver = new VNODE(ad);
  
  Solver->integrate(t,y,tend);
  
  if(!Solver->successful())
    {
      cerr << "VNODE-LP could not reach " 
	   << tend << endl;
      cerr << "Test 0-1 failed" << endl;
      return 1;
    }
  
  cout << "Test 0-1 SUCCESSFUL" << endl;
  return 0;
}

