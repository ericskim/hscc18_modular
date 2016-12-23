/*58:*/
#line 121 "./odectrl.w"

#include <fstream> 
#include "vnode.h"
using namespace std;
using namespace vnodelp;
/*51:*/
#line 9 "./odectrl.w"

struct
LorenzConsts{interval beta;double rho,sigma;};


/*:51*/
#line 126 "./odectrl.w"

/*53:*/
#line 27 "./odectrl.w"

template<typename var_type> 
void Lorenz2(int n,var_type*yp,const var_type*y,var_type t,
void*param){
LorenzConsts*p= (LorenzConsts*)param;
interval beta= p->beta;
double sigma= p->sigma;
double rho= p->rho;

yp[0]= sigma*(y[1]-y[0]);
yp[1]= y[0]*(rho-y[2])-y[1];
yp[2]= y[0]*y[1]-beta*y[2];
}





/*:53*/
#line 127 "./odectrl.w"


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
#line 131 "./odectrl.w"

/*52:*/
#line 15 "./odectrl.w"

LorenzConsts p;
p.sigma= 10.0;
p.beta= interval(8.0)/3.0;
p.rho= 28.0;



/*:52*/
#line 132 "./odectrl.w"

/*54:*/
#line 47 "./odectrl.w"

AD*ad= new FADBAD_AD(n,Lorenz2,Lorenz2,&p);

/*:54*/
#line 133 "./odectrl.w"

/*23:*/
#line 48 "./basic.w"

VNODE*Solver= new VNODE(ad);



/*:23*/
#line 134 "./odectrl.w"


/*55:*/
#line 58 "./odectrl.w"

Solver->setOneStep(on);
ofstream outFile1("odeparam1.out",ios::out);
while(t!=tend)
{
Solver->integrate(t,y,tend);
outFile1<<midpoint(y[0])<<"\t"
<<midpoint(y[1])<<"\t"
<<midpoint(y[2])<<endl;
}
outFile1.close();

/*:55*/
#line 136 "./odectrl.w"

/*57:*/
#line 103 "./odectrl.w"

/*56:*/
#line 82 "./odectrl.w"
t= 0.0;y[0]
= 15;y[1]= 15;y[2]= 36;

interval tend2= tend/2.0;

Solver->setFirstEntry();
while(t!=tend2)
Solver->integrate(t,y,tend2);


/*:56*/
#line 104 "./odectrl.w"


p.beta= 5.0;
ad->eval(&p);
ofstream outFile2("odeparam2.out",ios::out);
while(t!=tend)
{
Solver->integrate(t,y,tend);
outFile2<<midpoint(y[0])<<"\t"<<midpoint(y[1])<<"\t"
<<midpoint(y[2])<<endl;

}
outFile2.close();


/*:57*/
#line 137 "./odectrl.w"

return 0;
}


/*:58*/
