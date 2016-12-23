/*61:*/
#line 74 "./setintegp.w"

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
#line 81 "./setintegp.w"


int main()
{
/*21:*/
#line 24 "./basic.w"

const int n= 3;
interval t= 0.0,tend= 20.0;
iVector y(n);
y[0]= 15.0;
y[1]= 15.0;
y[2]= 36.0;


/*:21*/
#line 85 "./setintegp.w"

/*22:*/
#line 41 "./basic.w"

AD*ad= new FADBAD_AD(n,Lorenz,Lorenz);



/*:22*/
#line 86 "./setintegp.w"

/*23:*/
#line 48 "./basic.w"

VNODE*Solver= new VNODE(ad);



/*:23*/
#line 87 "./setintegp.w"

/*60:*/
#line 66 "./setintegp.w"

Solver->setTols(1e-14,1e-14);
Solver->setOrder(40);
Solver->setHmin(1e-5);


/*:60*/
#line 88 "./setintegp.w"

/*62:*/
#line 139 "./setintegp.w"

ofstream outFile1("lorenz.tight",ios::out);
ofstream outFile2("lorenz.step",ios::out);
ofstream outFile3("lorenz.apriori",ios::out);

/*:62*/
#line 89 "./setintegp.w"

/*43:*/
#line 58 "./basici.w"

Solver->setOneStep(on);

/*:43*/
#line 90 "./setintegp.w"

/*64:*/
#line 161 "./setintegp.w"

outFile1<<midpoint(t)<<"\t"
<<inf(y[0])<<"\t"<<sup(y[0])<<"\t"
<<width(y[0])<<endl;

/*:64*/
#line 91 "./setintegp.w"

while(t!=tend)
{
Solver->integrate(t,y,tend);
/*65:*/
#line 167 "./setintegp.w"


outFile1<<midpoint(t)<<"\t"
<<inf(y[0])<<"\t"<<sup(y[0])<<"\t"
<<width(y[0])<<endl;

outFile2<<midpoint(t)<<"\t"<<Solver->getStepsize()<<endl;


iVector Y= Solver->getAprioriEncl();
interval Tj= Solver->getT();

outFile3<<inf(Tj)<<"\t"<<inf(Y[0])<<endl;
outFile3<<inf(Tj)<<"\t"<<sup(Y[0])<<endl<<endl;

outFile3<<sup(Tj)<<"\t"<<inf(Y[0])<<endl;
outFile3<<sup(Tj)<<"\t"<<sup(Y[0])<<endl<<endl;

outFile3<<inf(Tj)<<"\t"<<inf(Y[0])<<endl;
outFile3<<sup(Tj)<<"\t"<<inf(Y[0])<<endl<<endl;

outFile3<<inf(Tj)<<"\t"<<sup(Y[0])<<endl;
outFile3<<sup(Tj)<<"\t"<<sup(Y[0])<<endl<<endl;




/*:65*/
#line 95 "./setintegp.w"

}
/*63:*/
#line 145 "./setintegp.w"

outFile1.close();
outFile2.close();
outFile3.close();








/*:63*/
#line 97 "./setintegp.w"

return 0;
}


/*:61*/
