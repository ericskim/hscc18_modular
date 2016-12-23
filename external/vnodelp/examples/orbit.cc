/*73:*/
#line 50 "./orbit.w"


#include <fstream> 
#include <sstream> 
#include <string> 
#include <cstdlib> 

#include "vnode.h"

using namespace std;
using namespace vnodelp;

template<typename var_type> 
void Orbit(int n,var_type*yp,const var_type*y,var_type t,
void*param)
{
interval mu= string_to_interval("0.012277471");
interval mu_h= 1.0-mu;

var_type D1= pow(sqr(y[0]+mu)+sqr(y[1]),interval(1.5));
var_type D2= pow(sqr(y[0]-mu_h)+sqr(y[1]),interval(1.5));

yp[0]= y[2];
yp[1]= y[3];

yp[2]= y[0]+2.0*y[3]-mu_h*(y[0]+mu)/D1-mu*(y[0]-mu_h)/D2;
yp[3]= y[1]-2.0*y[2]-mu_h*y[1]/D1-mu*y[1]/D2;
}


int main()
{
const int n= 4;
iVector y(n);

y[0]= string_to_interval("0.994");
y[1]= 0;
y[2]= 0;
y[3]= string_to_interval("-2.00158510637908252240537862224");
interval t= 0.0,tend= 35;

AD*ad= new FADBAD_AD(n,Orbit,Orbit);
VNODE*Solver= new VNODE(ad);

ofstream outFileSol("orbit_sol.out",ios::out);
ofstream outFileStep("orbit_step.out",ios::out);



outFileSol<<midpoint(y[0])<<"\t"<<midpoint(y[1])<<endl;
Solver->setOneStep(on);
while(t!=tend)
{
Solver->integrate(t,y,tend);
outFileSol<<midpoint(y[0])<<"\t"
<<midpoint(y[1])<<endl;
outFileStep<<midpoint(t)<<"\t"
<<Solver->getStepsize()<<endl;
}
outFileSol.close();
outFileStep.close();

return 0;
}




#line 1 "./vanderpol.w"

/*:73*/
