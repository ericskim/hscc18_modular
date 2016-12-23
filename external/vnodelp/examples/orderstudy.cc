/*70:*/
#line 95 "./orderstudy.w"

#include <fstream> 
#include <sstream> 
#include <string> 
#include <cstdlib> 

#include "vnode.h"
using namespace std;
using namespace vnodelp;
/*19:*/
#line 76 "./usage.w"

template<typename var_type> 
void Lorenz(int n,var_type*yp,
const var_type*y,
var_type t,void*param){

interval sigma(10.0),rho(28.0);
interval beta= interval(8.0)/3.0;

yp[0]= sigma*(y[1]-y[0]);
yp[1]= y[0]*(rho-y[2])-y[1];
yp[2]= y[0]*y[1]-beta*y[2];
}


#line 1 "./basic.w"

/*:19*/
#line 104 "./orderstudy.w"

/*68:*/
#line 34 "./orderstudy.w"

static double tol[]= {1e-7,1e-8,1e-9,1e-10,1e-11,1e-12,1e-13};
int main()
{
const int n= 3;
/*22:*/
#line 41 "./basic.w"

AD*ad= new FADBAD_AD(n,Lorenz,Lorenz);



/*:22*/
#line 39 "./orderstudy.w"

/*23:*/
#line 48 "./basic.w"

VNODE*Solver= new VNODE(ad);



/*:23*/
#line 40 "./orderstudy.w"


iVector y(n);

for(int i= 0;i<7;i++)
{
Solver->setTols(tol[i]);

/*69:*/
#line 84 "./orderstudy.w"

string prefix("order");


std::stringstream num(std::stringstream::out);
num<<tol[i];
string file_name= prefix+num.str()+".out";



/*:69*/
#line 48 "./orderstudy.w"

ofstream outFile(file_name.c_str(),ios::out);

cout<<" tol = "<<tol[i]<<": writing into "<<file_name<<" ..."<<endl;


for(int p= 5;p<=40;p++)

{
Solver->setOrder(p);
interval t= 0.0,tend= 10.0;

y[0]= 15.0;
y[1]= 15.0;
y[2]= 36.0;

Solver->setFirstEntry();
double time= getTime();

/*24:*/
#line 57 "./basic.w"

Solver->integrate(t,y,tend);


/*:24*/
#line 67 "./orderstudy.w"

/*25:*/
#line 62 "./basic.w"

if(!Solver->successful())
cout<<"VNODE-LP could not reach t = "<<tend<<endl;

/*:25*/
#line 68 "./orderstudy.w"

time= getTotalTime(time,getTime());
outFile<<p<<"\t"<<time<<endl;
}
outFile.close();
}

return 0;
}






/*:68*/
#line 105 "./orderstudy.w"



#line 1 "./sizestudy.w"

/*:70*/
