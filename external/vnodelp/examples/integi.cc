/*45:*/
#line 122 "./basici.w"

#include <fstream> 
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
#line 129 "./basici.w"


int main(){
/*42:*/
#line 14 "./basici.w"

/*21:*/
#line 24 "./basic.w"

const int n= 3;
interval t= 0.0,tend= 20.0;
iVector y(n);
y[0]= 15.0;
y[1]= 15.0;
y[2]= 36.0;


/*:21*/
#line 15 "./basici.w"

interval eps= interval(-1,1)/1e4;
y[0]+= eps;
y[1]+= eps;
y[2]+= eps;


/*:42*/
#line 132 "./basici.w"

/*22:*/
#line 41 "./basic.w"

AD*ad= new FADBAD_AD(n,Lorenz,Lorenz);



/*:22*/
#line 133 "./basici.w"

/*23:*/
#line 48 "./basic.w"

VNODE*Solver= new VNODE(ad);



/*:23*/
#line 134 "./basici.w"


/*44:*/
#line 70 "./basici.w"

/*43:*/
#line 58 "./basici.w"

Solver->setOneStep(on);

/*:43*/
#line 71 "./basici.w"

ofstream outFile("lorenzi.out",ios::out);


while(t!=tend)
{
Solver->integrate(t,y,tend);
if(Solver->successful()&&Solver->getGlobalExcess()<=15.0)
{
outFile<<midpoint(t)<<"\t"
<<midpoint(y[0])<<"\t"
<<rad(y[0])
<<"\t"<<
Solver->getGlobalExcess()<<endl;
}
else break;
}
outFile.close();







/*:44*/
#line 136 "./basici.w"

return 0;
}

/*:45*/
