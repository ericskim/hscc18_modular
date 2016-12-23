/*35:*/
#line 44 "./scalarp.w"

#include <ostream> 
#include "vnode.h"
using namespace std;
using namespace vnodelp;

/*32:*/
#line 15 "./scalarp.w"


template<typename var_type> 
void ScalarExample(int n,
var_type*yp,
const var_type*y,
var_type t,
void*param){

yp[0]= -y[0];
}


/*:32*/
#line 50 "./scalarp.w"


int main()
{
/*33:*/
#line 29 "./scalarp.w"


const int n= 1;
interval t= 0.0,tend= 20.0;
iVector y(n);
y[0]= 1.0;


/*:33*/
#line 54 "./scalarp.w"

/*34:*/
#line 38 "./scalarp.w"

AD*ad= new FADBAD_AD(n,ScalarExample,ScalarExample);


/*:34*/
#line 55 "./scalarp.w"

/*23:*/
#line 48 "./basic.w"

VNODE*Solver= new VNODE(ad);



/*:23*/
#line 56 "./scalarp.w"

/*24:*/
#line 57 "./basic.w"

Solver->integrate(t,y,tend);


/*:24*/
#line 57 "./scalarp.w"

/*25:*/
#line 62 "./basic.w"

if(!Solver->successful())
cout<<"VNODE-LP could not reach t = "<<tend<<endl;

/*:25*/
#line 58 "./scalarp.w"

/*26:*/
#line 67 "./basic.w"

cout<<"Solution enclosure at t = "<<t<<endl;
printVector(y);




/*:26*/
#line 59 "./scalarp.w"

return 0;
}


/*:35*/
