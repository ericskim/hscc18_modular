/*27:*/
#line 80 "./basic.w"

#include <ostream> 
#include "vnode.h"
using namespace std;
using namespace vnodelp;

/*20:*/
#line 6 "./basic.w"

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
#line 7 "./basic.w"

int main(){

/*21:*/
#line 24 "./basic.w"

const int n= 3;
interval t= 0.0,tend= 20.0;
iVector y(n);
y[0]= 15.0;
y[1]= 15.0;
y[2]= 36.0;


/*:21*/
#line 10 "./basic.w"

/*22:*/
#line 41 "./basic.w"

AD*ad= new FADBAD_AD(n,Lorenz,Lorenz);



/*:22*/
#line 11 "./basic.w"

/*23:*/
#line 48 "./basic.w"

VNODE*Solver= new VNODE(ad);



/*:23*/
#line 12 "./basic.w"

/*24:*/
#line 57 "./basic.w"

Solver->integrate(t,y,tend);


/*:24*/
#line 13 "./basic.w"

/*25:*/
#line 62 "./basic.w"

if(!Solver->successful())
cout<<"VNODE-LP could not reach t = "<<tend<<endl;

/*:25*/
#line 14 "./basic.w"

/*26:*/
#line 67 "./basic.w"

cout<<"Solution enclosure at t = "<<t<<endl;
printVector(y);




/*:26*/
#line 15 "./basic.w"

return 0;
}

/*:20*/
#line 86 "./basic.w"




/*:27*/
