#include <ostream> 
#include "vnode.h"

using namespace std;
using namespace vnodelp;

template<typename var_type> 
void P6(int n, var_type*yp, const var_type*y, var_type t,
	void*param)
{
  var_type t2 = sqr(t);
  
  yp[0] = sin(t+10.0)*y[0] -  2.0*y[1] - y[2]  -cos(t);
  
  yp[1] = 3.0*y[0] - 4.0*cos(t2)*y[1]          -cos(t);
  
  var_type et2 = exp(-t2);
  yp[2] = et2*y[0] - et2*y[1]                  -sin(t);
}

int main()
{
  const int n = 3;
  
  interval t = 0.0, tend = 20;
  iVector y(n);
  
  y[0] = interval(0,5);
  y[1] = interval(-2,6);
  y[2] = interval(5,12);

  AD *ad = new FADBAD_AD(n,P6,P6);
  VNODE *Solver = new VNODE(ad);
  
  // integrate with the above initial condition and save the result.
  Solver->integrate(t,y,tend);
  iVector yi = y;
  
  // integrate with each corner point
  // 1
  iVector maple_solution(3);
  maple_solution[0] = string_to_interval("73.425783284967972375");
  maple_solution[1] = string_to_interval("-35.161996476428791191");
  maple_solution[2] = string_to_interval("5.0809398723127131529");

  t = 0; tend = 20;
  y[0] = 0; y[1] = -2; y[2] = 5;
  Solver->setFirstEntry();
  Solver->integrate(t,y,tend);
  if (disjoint(y,yi) || !subseteq(maple_solution,y ))
    {
      cerr << "Test 6 failed" << endl;
      return 1;
    }
  
  // 2
  maple_solution[0] = string_to_interval("163.80316188252839853");
  maple_solution[1] = string_to_interval("-77.987614188132603485");
  maple_solution[2] = string_to_interval("12.057253902743765676");

  y[0] = 0; y[1] = -2; y[2] = 12;
  t = 0; tend = 20;
  Solver->setFirstEntry();
  Solver->integrate(t,y,tend);
  if ( disjoint(y,yi) || !subseteq(maple_solution,y) )
    {
      cerr << "Test 6 failed" << endl;
      return 1;
    }

  // 3
  maple_solution[0] = string_to_interval("55.361087355369320947");
  maple_solution[1] = string_to_interval("-26.014987809924523279");
  maple_solution[2] = string_to_interval("2.8533625088266014731");

  y[0] = 0; y[1] = 6; y[2] = 5;
  t = 0; tend = 20;
  Solver->setFirstEntry();
  Solver->integrate(t,y,tend);
  if (disjoint(y,yi)  || !subseteq(maple_solution,y) )
    {
      cerr << "Test 6 failed" << endl;
      return 1;
    }

  // 4 
  maple_solution[0] = string_to_interval("145.73846595292974710");
  maple_solution[1] = string_to_interval("-68.840605521628335573");
  maple_solution[2] = string_to_interval("9.8296765392576539965");

  y[0] = 0.0; y[1] = 6.0; y[2] = 12.0;
  t = 0; tend = 20;
  Solver->setFirstEntry();
  Solver->integrate(t,y,tend);
  if (disjoint(y,yi) )
    {
      cerr << "Test 6 failed" << endl;
      return 1;
    }
  
  
  // 5
  maple_solution[0] = string_to_interval("80.110231016749208330");
  maple_solution[1] = string_to_interval("-38.548283829601594817");
  maple_solution[2] = string_to_interval("5.7339951957821465960");

  y[0] = 5; y[1] = -2; y[2] = 5;
  t = 0; tend = 20;
  Solver->setFirstEntry();
  Solver->integrate(t,y,tend);
  if (disjoint(y,yi)  || !subseteq(maple_solution,y) )
    {
      cerr << "Test 6 failed" << endl;
      return 1;
    }
  
  // 6 
  maple_solution[0] = string_to_interval("170.48760961430963449");
  maple_solution[1] = string_to_interval("-81.373901541305407111");
  maple_solution[2] = string_to_interval("12.710309226213199119");

  y[0] = 5; y[1] = -2; y[2] = 12;
  t = 0; tend = 20;
  Solver->setFirstEntry();
  Solver->integrate(t,y,tend);
  if ( disjoint(y,yi) || !subseteq(maple_solution,y) ) 
    {
      cerr << "Test 6 failed" << endl;
      return 1;
    }
  
  // 7
  maple_solution[0] = string_to_interval("62.045535087150556901");
  maple_solution[1] = string_to_interval("-29.401275163097326905");
  maple_solution[2] = string_to_interval("3.5064178322960349163");

  y[0] = 5; y[1] = 6; y[2] = 5;
  t = 0; tend = 20;
  Solver->setFirstEntry();
  Solver->integrate(t,y,tend);
  if (disjoint(y,yi))
    {
      cerr << "Test 6 failed" << endl;
      return 1;
    }
  
  // 8 
  maple_solution[0] = string_to_interval("152.42291368471098306");
  maple_solution[1] = string_to_interval("-72.226892874801139199");
  maple_solution[2] = string_to_interval("10.482731862727087440");

  y[0] = 5; y[1] = 6; y[2] = 12;
  t = 0; tend = 20;
  Solver->setFirstEntry();
  Solver->integrate(t,y,tend);
  if ( disjoint(y,yi) || !subseteq(maple_solution,y) ) 
    {
      cerr << "Test 6 failed" << endl;
      return 1;
    }
  
  cout <<"Test 6 SUCCESSFUL\n";
  return 0;
}

