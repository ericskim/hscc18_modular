#include <ostream> 
#include "vnode.h"

using namespace std;
using namespace vnodelp;

#define N_SIZE 3
  

#define X(i)  Y[i-1]
#define Y(i)  Y[N_SIZE+i-1]
#define VX(i) Y[2*N_SIZE+i-1]
#define VY(i) Y[3*N_SIZE+i-1]

#define XP(i)  Yp[i-1]
#define YP(i)  Yp[N_SIZE+i-1]
#define VXP(i) Yp[2*N_SIZE+i-1]
#define VYP(i) Yp[3*N_SIZE+i-1]


template<typename var_type> 
void ThreeBody(int n, var_type*Yp, const var_type*Y, var_type t,
	       void*param)
{
  for ( int i = 1; i <= N_SIZE; i++ )
    {
      XP(i) = VX(i);
      YP(i) = VY(i);
      
      VXP(i) = interval(0.0);
      VYP(i) = interval(0.0);
      
      for ( int j = 1; j <= N_SIZE; j++ )
	if ( j != i )
	  {
	    var_type tmp = sqr( X(j) - X(i) ) + sqr( Y(j) - Y(i) );
	    tmp *= sqrt(tmp);
	    
	    VXP(i) += (X(j)-X(i))/tmp;
	    VYP(i) += (Y(j)-Y(i))/tmp;
	  }
    }
}



int main()
{
  const int n = 12;
  iVector y(n);
  
  // x
  y[0] = 0.0;
  y[1] = 0.0;
  y[2] = 0.0;
  
  // y 
  y[3] = -0.070946548000287;
  y[4] =  1.075903540029147;
  y[5] = -1.004956992028862;
  
  // dx
  y[6] =  1.231870762148251;
  y[7] = -0.195089035291636;
  y[8] = -1.036781726856616;
  
  // dy
  y[9]  = 0.0;
  y[10] = 0.0;
  y[11] = 0.0;
  
  iVector y0 = y; // save IC

  AD *ad = new FADBAD_AD(n,ThreeBody,ThreeBody);
  VNODE *Solver = new VNODE(ad);
  
  Solver->setTols(1e-15,1e-15);
  
  // integrate from 0 to 1
  interval t = 0.0, tend = 1.0;
  Solver->integrate(t,y,tend);

  // integrate from 1 to 0
  tend = 0;
  Solver->setFirstEntry();
  Solver->integrate(t,y,tend);
  
  if (!subseteq(y0,y))
    {
      cerr << "Test n-2 failed: y0 must be contained in y" << endl;
      printVector(y0, "y0 ");
      printVector(y,  "y  ");
      return 1;
    }
  
  cout << "Test n-2 SUCCESSFUL" << endl;
  return 0;
}

