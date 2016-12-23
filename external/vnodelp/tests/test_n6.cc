#include <ostream> 
#include "vnode.h"

using namespace std;
using namespace vnodelp;

#define y(i)   Y[i-1]
#define f(i)  Yp[i-1]

template<typename var_type> 
void Pleiades(int n,var_type*Yp, const var_type*Y, var_type t,
	      void*param)
{
  for (int i = 1; i<=7; i++)
    {
      var_type sumx = interval(0.0);
      var_type sumy = interval(0.0);
      
      for (int j = 1; j<=7; j++)
	{
	  int mj = j;
	  var_type rij = sqr( y(i)-y(j) ) + sqr( y(i+7)-y(j+7));
	  var_type rij32 = rij*sqrt(rij);
	  if (j != i) 
	    {
	      sumx = sumx+(double)mj*(y(j)-y(i))/rij32;
	      sumy = sumy+(double)mj*(y(j+7)-y(i+7))/rij32;
	    }
	}
      
      f(i+14) = sumx;
      f(i+21) = sumy;
    }
  
  for (int i = 1; i <= 14; i++ )
    f(i) = y(i+14);
}


int main()
{
  const int n = 28;
  iVector Y(n);
  
  y(1)  =  3.0;
  y(2)  =  3.0;
  y(3)  = -1.0;
  y(4)  = -3.0;
  y(5)  =  2.0;
  y(6)  = -2.0;
  y(7)  =  2.0;
  y(8)  =  3.0;
  y(9)  = -3.0;
  y(10) =  2.0;
  y(11) =  0.0;
  y(12) =  0.0;
  y(13) = -4.0;
  y(14) =  4.0;
  y(15) =  0.0;
  y(16) =  0.0;
  y(17) =  0.0;
  y(18) =  0.0;
  y(19) =  0.0;
  y(20) =  1.75;
  y(21) = -1.5;
  y(22) =  0.0;
  y(23) =  0.0;
  y(24) =  0.0;
  y(25) = -1.25;
  y(26) =  1.0;
  y(27) =  0.0;
  y(28) =  0.0;
  
  interval t = 0.0, tend = 3;

  AD *ad = new FADBAD_AD(n,Pleiades,Pleiades);
  VNODE *Solver = new VNODE(ad);
  
  Solver->setTols(1e-16,1e-16);
  Solver->setOrder(15);
  
  Solver->setOneStep(on);
  int step_count = 0;
  while ( t!=tend )
    {
      step_count++;
      Solver->integrate(t,Y,tend);
      if ( (step_count % 10 ) == 0 )
	cout << "Step number " << step_count 
	     << " Computed at t = " << t 
	     << "  Excess = " << Solver->getGlobalExcess() << endl;
    }

  
  iVector ref_sol(n);
  
#define yr(i)   ref_sol[i-1]
  
  yr(1) =  0.3706139143970502E+000;
  yr(2) =  0.3237284092057233E+001;
  yr(3) = -0.3222559032418324E+001;
  yr(4) =  0.6597091455775310E+000;
  yr(5) =  0.3425581707156584E+000;
  yr(6) =  0.1562172101400631E+001;
  yr(7) = -0.7003092922212495E+000;
  yr(8) = -0.3943437585517392E+001;
  yr(9) = -0.3271380973972550E+001;
  yr(10) =  0.5225081843456543E+001;
  yr(11) = -0.2590612434977470E+001;
  yr(12) =  0.1198213693392275E+001;
  yr(13) = -0.2429682344935824E+000;
  yr(14) =  0.1091449240428980E+001;
  yr(15) =  0.3417003806314313E+001;
  yr(16) =  0.1354584501625501E+001;
  yr(17) = -0.2590065597810775E+001;
  yr(18) =  0.2025053734714242E+001;
  yr(19) = -0.1155815100160448E+001;
  yr(20) = -0.8072988170223021E+000;
  yr(21) =  0.5952396354208710E+000;
  yr(22) = -0.3741244961234010E+001;
  yr(23) =  0.3773459685750630E+000;
  yr(24) =  0.9386858869551073E+000;
  yr(25) =  0.3667922227200571E+000;
  yr(26) = -0.3474046353808490E+000;
  yr(27) =  0.2344915448180937E+001;
  yr(28) = -0.1947020434263292E+001;
  
  
  /* We subtract the reference solution from the IVP Test Set from the

     computed bounds. */
  subViVi(Y,ref_sol);
  if ( inf_normV(Y) > 1e-2)
    {
      cerr << "Test n-6 FAILED" << endl;
      return 1;
    }
  
  cout << "Test n-6 SUCCESSFUL" << endl;
  
  return 0;
}




