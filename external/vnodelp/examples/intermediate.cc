/*48:*/
#line 30 "./intermed.w"

#include <iostream> 
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
#line 35 "./intermed.w"


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
#line 38 "./intermed.w"

/*22:*/
#line 41 "./basic.w"

AD*ad= new FADBAD_AD(n,Lorenz,Lorenz);



/*:22*/
#line 39 "./intermed.w"

/*23:*/
#line 48 "./basic.w"

VNODE*Solver= new VNODE(ad);



/*:23*/
#line 40 "./intermed.w"

/*47:*/
#line 9 "./intermed.w"


interval step= string_to_interval("0.1");


tend= 0.0;
for(int i= 1;i<=3;i++)
{
tend+= step;
Solver->integrate(t,y,tend);
/*26:*/
#line 67 "./basic.w"

cout<<"Solution enclosure at t = "<<t<<endl;
printVector(y);




/*:26*/
#line 19 "./intermed.w"

}

tend= 10;
Solver->integrate(t,y,tend);
/*26:*/
#line 67 "./basic.w"

cout<<"Solution enclosure at t = "<<t<<endl;
printVector(y);




/*:26*/
#line 24 "./intermed.w"




/*:47*/
#line 41 "./intermed.w"

return 0;
}

/*:48*/
