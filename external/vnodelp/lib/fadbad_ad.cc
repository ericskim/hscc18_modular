/*350:*/
#line 169 "./fadbad.w"

#include "fadbad_ad.h"
namespace vnodelp{

/*340:*/
#line 52 "./fadbad.w"

FadbadODE::FadbadODE(int n,Tfunction f,void*param):AD_ODE()
{
size= n;
y_coeff= new Tinterval[2*n];
f_coeff= y_coeff+n;
fcn= f;
fcn(size,f_coeff,y_coeff,t,param);
}

/*:340*/
#line 173 "./fadbad.w"

/*341:*/
#line 63 "./fadbad.w"

FadbadODE::~FadbadODE()
{
delete y_coeff;
}



/*:341*/
#line 174 "./fadbad.w"

/*342:*/
#line 72 "./fadbad.w"

void FadbadODE::set(const interval&t0,const iVector&y0,
const interval&h,int k)
{
t[0]= t0;
t[1]= h;
stepsize= h;order= k;
for(int eqn= 0;eqn<size;eqn++)
y_coeff[eqn][0]= y0[eqn];
}



/*:342*/
#line 175 "./fadbad.w"

/*343:*/
#line 90 "./fadbad.w"

void FadbadODE::compTerms()
{


for(int eqn= 0;eqn<size;eqn++)
f_coeff[eqn].reset();


for(int coeff= 1;coeff<=order;coeff++)

for(int eqn= 0;eqn<size;eqn++)

{
f_coeff[eqn].eval(coeff-1);
y_coeff[eqn][coeff]= stepsize*f_coeff[eqn][coeff-1]/double(coeff);

}

}


/*:343*/
#line 176 "./fadbad.w"

/*344:*/
#line 112 "./fadbad.w"

void FadbadODE::sumTerms(iVector&sum,int m)
{
interval s;
for(int eqn= 0;eqn<size;eqn++)
{
s= 0.0;
for(int coeff= m;coeff>=0;coeff--)
s+= y_coeff[eqn][coeff];
sum[eqn]= s;
}
}


/*:344*/
#line 177 "./fadbad.w"

/*345:*/
#line 126 "./fadbad.w"

void FadbadODE::getTerm(iVector&term,int i)const
{
for(int eqn= 0;eqn<size;eqn++)
term[eqn]= y_coeff[eqn][i];
}

/*:345*/
#line 178 "./fadbad.w"

/*346:*/
#line 134 "./fadbad.w"

interval FadbadODE::getStepsize()const{return stepsize;}



/*:346*/
#line 179 "./fadbad.w"

/*347:*/
#line 140 "./fadbad.w"

void FadbadODE::eval(void*param){fcn(size,f_coeff,y_coeff,t,param);}




/*:347*/
#line 180 "./fadbad.w"

}







/*:350*/
