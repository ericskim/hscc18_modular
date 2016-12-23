/*71:*/
#line 22 "./sizestudy.w"

#include <ostream> 
#include "vnode.h"

using namespace std;
using namespace vnodelp;

template<typename var_type> 
void DETEST_C3(int n,var_type*yp,const var_type*y,var_type t,
void*param)
{
yp[0]= -2.0*y[0]+y[1];

for(int i= 1;i<n-1;i++)
{
yp[i]= y[i-1]-2.0*y[i]+y[i+1];
}

yp[n-1]= y[n-2]-2.0*y[n-1];
}


int main()
{
for(int n= 40;n<=200;n+= 20)
{
cout<<n;

interval t= 0.0,tend= 5;
iVector y(n);

for(int i= 0;i<n;i++)
y[i]= 0;
y[0]= 1;

AD*ad= new FADBAD_AD(n,DETEST_C3,DETEST_C3);
VNODE*Solver= new VNODE(ad);

double time_start= getTime();

Solver->integrate(t,y,tend);

double time_end= getTime();
cout<<"  "<<getTotalTime(time_start,time_end)
<<"   "<<Solver->steps
<<endl;

delete Solver;
delete ad;
}
return 0;
}


/*:71*/
