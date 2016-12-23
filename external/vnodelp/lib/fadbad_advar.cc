/*361:*/
#line 358 "./fadbad.w"

#include "fadbad_advar.h"
namespace vnodelp{
/*353:*/
#line 232 "./fadbad.w"

FadbadVarODE::FadbadVarODE(int n,TFfunction f,void*param)


{
size= n;
tf_in= new TFinterval[2*n];
tf_out= tf_in+n;
fcn= f;
fcn(size,tf_out,tf_in,tf,param);
}

/*:353*/
#line 361 "./fadbad.w"

/*354:*/
#line 245 "./fadbad.w"

FadbadVarODE::~FadbadVarODE()
{
delete[]tf_in;
}

/*:354*/
#line 362 "./fadbad.w"

/*355:*/
#line 252 "./fadbad.w"

void FadbadVarODE::set(
const interval&t0,const iVector&y0,const interval&h,int k)
{
stepsize= h;
order= k;
tf[0].x()= t0;
tf[1].x()= h;

for(int eqn= 0;eqn<size;eqn++)
{
tf_in[eqn][0]= y0[eqn];
}
}



/*:355*/
#line 363 "./fadbad.w"

/*356:*/
#line 270 "./fadbad.w"

void FadbadVarODE::compTerms()
{
for(int eqn= 0;eqn<size;eqn++)
tf_out[eqn].reset();

for(int eqn= 0;eqn<size;eqn++)
tf_in[eqn][0].diff(eqn,size);


for(int coeff= 0;coeff<order;coeff++)
{
for(int eqn= 0;eqn<size;eqn++)
{
tf_out[eqn].eval(coeff);

tf_in[eqn][coeff+1]= stepsize*(tf_out[eqn][coeff]/double(coeff+1));
}
}
}


/*:356*/
#line 364 "./fadbad.w"

/*357:*/
#line 293 "./fadbad.w"



void FadbadVarODE::sumTerms(iMatrix&Sum,int k)
{

for(int row= 0;row<size;row++)
for(int col= 0;col<size;col++)
{
interval s= 0.0;

for(int coeff= k;coeff>=1;coeff--)s+= tf_in[row][coeff].d(col);
Sum[row][col]= s;

}

for(int row= 0;row<size;row++)
Sum[row][row]+= 1.0;
}


/*:357*/
#line 365 "./fadbad.w"

/*358:*/
#line 314 "./fadbad.w"

void FadbadVarODE::getTerm(iMatrix&Term,int i)const
{
for(int row= 0;row<size;row++)
for(int col= 0;col<size;col++)
Term[row][col]= tf_in[row][i].d(col);
}










/*:358*/
#line 366 "./fadbad.w"

}

/*:361*/
