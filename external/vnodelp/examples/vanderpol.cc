/*74:*/
#line 54 "./vanderpol.w"

#include <fstream> 
#include <iomanip> 
#include <sstream> 
#include <string> 


#include "vnode.h"

using namespace std;
using namespace vnodelp;

template<typename var_type> 
void VDP(int n,var_type*yp,const var_type*y,var_type t,
void*param)
{
double*MU= (double*)param;
yp[0]= y[1];
 yp[1]= (*MU)*(1.0-sqr(y[0]))*y[1]-y[0];
}




int main()
{
const int n= 2;

double MU= 10.0;
AD*ad= new FADBAD_AD(n,VDP,VDP,&MU);
VNODE*Solver= new VNODE(ad);


ofstream outSteps("vdp_nosteps.out",ios::out);
outSteps<<fixed<<showpoint<<setprecision(1);

for(int i= 1;i<=4;i++)
{


string prefix("vdp_step");
std::stringstream num(std::stringstream::out);
num<<i;
string file_name= prefix+num.str()+".out";
ofstream outStepSizes(file_name.c_str(),ios::out);

cout<<" MU = "<<MU<<" writing into "
<<file_name<<"..."<<endl;

Solver->setFirstEntry();
Solver->setOneStep(on);

interval t= 0,tend= 200;
iVector y(n);
y[0]= 2.0,y[1]= 0;

interval t_prev= t;
double time= getTime();
while(t!=tend)
{
Solver->integrate(t,y,tend);
if(midpoint(t-t_prev)>=0.01||t==tend)
{
outStepSizes<<midpoint(t)<<"\t"
<<Solver->getStepsize()<<endl;
t_prev= t;
}
}
outStepSizes.close();
time= getTotalTime(time,getTime());

outSteps<<"$10^{"<<i<<"}$"<<"\t & "
<<Solver->getNoSteps()
<<"\t & "
<<time<<"\t"
<<"\\\\"<<endl;

MU*= 10.0;
ad->eval(&MU);
}
outSteps.close();

return 0;
}

#line 79 "./vnode.w"
/*:74*/
