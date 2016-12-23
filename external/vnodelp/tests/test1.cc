#include <ostream> 
#include "vnode.h"

using namespace std;
using namespace vnodelp;

template<typename var_type> 
void Oscill(int n, var_type*yp, const var_type*y, var_type t,
	    void*param)
{
  yp[0] =  y[1];
  yp[1] = -y[0];
}


static void trueSolution( iVector &y, 
			  const interval & t0, const iVector &y0)
{
  y[0] =  cos(t0)*y0[0] + sin(t0)*y0[1];
  y[1] = -sin(t0)*y0[0] + cos(t0)*y0[1];
}


int main()
{
  const int n = 2;
  
  interval t = 0.0, tend = 10000;
  iVector y(n);
  y[0] = 1.0;
  y[1] = 1.0;
  
  AD *ad = new FADBAD_AD(n,Oscill,Oscill);
  VNODE *Solver = new VNODE(ad);
  
  iVector  y0 = y;
  iVector yt(n), tmp(n);  
  
  Solver->setOneStep(on);

  while ( t!=tend )
    {
      Solver->integrate(t,y,tend);
      if(!Solver->successful())
	{
	  cerr << "VNODE-LP could not reach " 
	       << tend << endl;
	  return 1;
	}
      
      trueSolution(yt,t,y0);
      
      if (!intersect(tmp,yt,y))
	{
	  cerr << "Test 1 failed: y and yt must intersect" << endl;
	  printVector(y, "y  ");
	  printVector(yt,"yt ");
	  return 1;
	}
    }
  
  cout << "Test 1 SUCCESSFUL" << endl;

  return 0;
}

