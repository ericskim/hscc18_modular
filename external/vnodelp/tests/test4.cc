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


int main()
{
  const int n = 2;
  interval t = 0.0;
  
  iVector y(n);
  y[0] = 1.0;
  y[1] = 1.0;
  iVector y0 = y;
  
  AD *ad = new FADBAD_AD(n, Oscill,Oscill);
  VNODE *Solver = new VNODE(ad);
  
  for ( int i = 2; i <= 1000; i+=2 )
    {
      interval tend = -double(i)*pi();
      Solver->integrate(t,y,tend);
      
      if (!subseteq(y0,y))
	{
	  cerr << "Test 4 failed: y0 must be contained in y" << endl;
	  printVector(y0, "y0 ");
	  printVector(y,  "y  ");
	  return 1;
	}
    }

  cout << "Test 4 SUCCESSFUL" << endl;
  return 0;
}

