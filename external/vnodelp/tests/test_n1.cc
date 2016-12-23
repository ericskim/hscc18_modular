#include <ostream> 
#include "vnode.h"

using namespace std;
using namespace vnodelp;

template<typename var_type> 
void Lorenz(int n, var_type*yp, const var_type*y, var_type t,
	    void*param)
{
  interval sigma(10.0),rho(28.0);
  interval beta = interval(8.0)/3.0;
  
  yp[0] = sigma*(y[1]-y[0]);
  yp[1] = y[0]*(rho-y[2])-y[1];
  yp[2] = y[0]*y[1]-beta*y[2];
}


int main()
{
  const int n = 3;
  
  iVector y(n);
  y[0] = 15.0;
  y[1] = 15.0;
  y[2] = 36.0;
  
  // save the initial condition
  iVector y0 = y;
  
  AD *ad = new FADBAD_AD(n,Lorenz,Lorenz);
  VNODE *Solver = new VNODE(ad);
  
  Solver->setTols(1e-15,1e-15);
  
  // integrate from 0 to 1
  interval t = 0.0, tend = 1.0;
  Solver->integrate(t,y,tend);
  
  // check if approximate solution is contained in the computed bounds
  iVector maple_solution(n);
  // approximate solution at t = 1
  maple_solution[0] = string_to_interval("-6.9453541599034593196");
  maple_solution[1] = string_to_interval("2.9971546266290307395");
  maple_solution[2] = string_to_interval("35.144350305722419178");
  
  if (!subseteq(maple_solution,y))
    {
      cerr << "Test n-1 failed: approximate solution must be contained in y" 
	   << endl;
      printVector(y);
      printVector(maple_solution);
      return 1;
    }
  
  // enclosure at 1
  iVector y1 = y; 
  
  assert(t==tend);
  tend = 0.0;

  // integrate from 1 to 0
  Solver->setFirstEntry();
  Solver->integrate(t,y,tend);
  // check if initial condition is in y0
  if (!subseteq(y0,y))
    {
      cerr << "Test n-1 failed: y0 must be contained in y" << endl;
      return 1;
    }
  
  // integrate from 0 to 1
  tend = 1.0;
  Solver->setFirstEntry();
  Solver->integrate(t,y,tend);
  
  // check if approximate solution is contained in y
  if (!subseteq(maple_solution,y))
    {
      cerr << "Test n-1 failed: approximate solution must be contained in y" 
	   << endl;
      return 1;
    }

  
  // check if previous and current enclosures at t = 1 intersect
  iVector tmp(n);
  if (!intersect(tmp, y1,y) )
    {
      cerr << "Test n-1 failed: y and y1 must intersect" << endl;
      return 1;
    }
  
  cout << "Test n-1 SUCCESSFUL" << endl;
  return 0;
}

