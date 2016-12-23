#include <ostream> 
#include "vnode.h"

using namespace std;
using namespace vnodelp;

template<typename var_type> 
void ExpDecay(int n, var_type*yp, const var_type*y, var_type t,
	    void*param)
{
  yp[0] =      y[0] - 2.0*y[1];
  yp[1] =  3.0*y[0] - 4.0*y[1];
}

void static trueSolution( iVector &y, 
			  const interval & t0, const iVector &y0)
{
  interval t1 = 5.0*exp(-t0);
  interval t2 = 4.0*exp(-2.0*t0);
  interval t3 = 6.0*exp(-2.0*t0);

  y[0] =  t1 - t2;
  y[1] =  t1 - t3;
}


int main()
{
  const int n = 2;
  
  interval t = 0.0, tend = 10000;
  iVector y(n);
  y[0] = 1.0;
  y[1] = -1.0;
  
  AD *ad = new FADBAD_AD(n,ExpDecay,ExpDecay);
  VNODE *Solver = new VNODE(ad);
  
  iVector  y0 = y;
  iVector yt(n), tmp(n);  
  
  Solver->setOneStep(on);
  
  while ( t!=tend )
    {
      Solver->integrate(t,y,tend);
      if(!Solver->successful())
	{
	  cerr << "Test 5 failed: " << endl;
	  cerr << "VNODE-LP could not reach " 
	       << tend << endl;
	  return 1;
	}
  
      trueSolution(yt,t,y0);
      if (!intersect(tmp,yt,y))
	{
	  cerr << "Test 5 failed: y and yt must intersect" << endl;
	  printVector(y, "y  ");
	  printVector(yt,"yt ");
	  return 1;
	}
    }
  cout << "Test 5 SUCCESSFUL" << endl;
  return 0;
}

