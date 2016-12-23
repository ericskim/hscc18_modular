/*40:*/
#line 58 "./timedep.w"

#include <ostream> 
#include "vnode.h"
using namespace std;
using namespace vnodelp;
/*37:*/
#line 20 "./timedep.w"


template<typename var_type> 
void DETEST_E1(int n,
var_type*yp,const var_type*y,
var_type t,void*param)
{

var_type t1= t+1.0;
yp[0]= y[1];
yp[1]= -(y[1]/t1+(1.0-0.25/(t1*t1))*y[0]);
}

/*:37*/
#line 63 "./timedep.w"


int main(){
/*38:*/
#line 38 "./timedep.w"


const int n= 2;
interval t= 0.0,tend= 20.0;
iVector y(n);
y[0]= string_to_interval("0.6713967071418030");
y[1]= string_to_interval("0.09540051444747446");



/*:38*/
#line 66 "./timedep.w"

/*39:*/
#line 49 "./timedep.w"

AD*ad= new FADBAD_AD(n,DETEST_E1,DETEST_E1);




/*:39*/
#line 67 "./timedep.w"

/*23:*/
#line 48 "./basic.w"

VNODE*Solver= new VNODE(ad);



/*:23*/
#line 68 "./timedep.w"

/*24:*/
#line 57 "./basic.w"

Solver->integrate(t,y,tend);


/*:24*/
#line 69 "./timedep.w"

/*25:*/
#line 62 "./basic.w"

if(!Solver->successful())
cout<<"VNODE-LP could not reach t = "<<tend<<endl;

/*:25*/
#line 70 "./timedep.w"

/*26:*/
#line 67 "./basic.w"

cout<<"Solution enclosure at t = "<<t<<endl;
printVector(y);




/*:26*/
#line 71 "./timedep.w"

return 0;
}


/*:40*/
