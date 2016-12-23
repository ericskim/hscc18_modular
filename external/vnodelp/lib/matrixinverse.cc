/*173:*/
#line 422 "./matrixinverse.w"

#include <limits> 

#include  "matrixinverse.h"
#include "basiclinalg.h"
#include "intvfuncs.h"
#include "debug.h"

using namespace std;

namespace v_blas{
/*171:*/
#line 377 "./matrixinverse.w"

MatrixInverse::MatrixInverse(int n)
{
M= new double[n*n];
ipiv= new int[n];
lwork= 2*n;
work= new double[lwork];

sizeV(radx,n);
sizeV(b0,n);
sizeV(x1,n);
sizeV(x,n);

sizeM(B,n);
sizeM(Ci,n);
sizeM(C,n);
}

MatrixInverse::~MatrixInverse(){
delete[]work;
delete[]ipiv;
delete[]M;
}


/*:171*/
#line 432 "./matrixinverse.w"

/*158:*/
#line 71 "./matrixinverse.w"

extern"C"
{

void dgetrf_(int*m,int*n,double*A,int*lda,
int*ipiv,int*info);

void dgetri_(int*n,double*A,int*lda,int*ipiv,double*work,
int*lwork,int*info);
}



bool MatrixInverse::invertMatrix(pMatrix&Ainv,const pMatrix&A)
{
int n= sizeM(A);
int lda= n;
int info;

matrix2pointer(M,A);











v_bias::round_nearest();

dgetrf_(&n,&n,M,&lda,ipiv,&info);

if(info!=0)
{
#ifdef VNODE_DEBUG
printMessage("Could not invert a matrix");
#endif
return false;
}









dgetri_(&n,M,&lda,ipiv,work,&lwork,&info);
if(info!=0)
{
#ifdef VNODE_DEBUG
printMessage("Could not invert a matrix");
#endif
return false;
}
pointer2matrix(Ainv,M);

return true;
}




/*:158*/
#line 433 "./matrixinverse.w"

/*160:*/
#line 175 "./matrixinverse.w"

bool
MatrixInverse::encloseLS(iVector&x,const pMatrix&A,const
iVector&b0,const iMatrix&B,double beta){








/*161:*/
#line 201 "./matrixinverse.w"

v_bias::round_down();
double a= 1-beta;
double a1= inf_normV(b0);
v_bias::round_up();
a= a1/a;


setV(x,v_bias::interval(-a,a));


/*:161*/
#line 187 "./matrixinverse.w"

/*165:*/
#line 248 "./matrixinverse.w"

double sum_old_radii= numeric_limits<double> ::max();
double sum_radii;

/*164:*/
#line 237 "./matrixinverse.w"

rad(radx,x);
v_bias::round_nearest();
sum_radii= accumulate(radx.begin(),radx.end(),0.0);


/*:164*/
#line 252 "./matrixinverse.w"

round_up();
int max_iterations= 20;
int counter= 0;
double mult= (1+beta)/2;

while(sum_radii<mult*sum_old_radii&&counter<max_iterations)
{
/*162:*/
#line 216 "./matrixinverse.w"

multMiVi(x1,B,x);
addViVi(x1,b0);

/*:162*//*163:*/
#line 222 "./matrixinverse.w"

bool b= intersect(x,x,x1);
if(b==false)
{
#ifdef VNODE_DEBUG
printMessage("x and x1 do not intersect");
#endif
return false;
}



/*:163*/
#line 260 "./matrixinverse.w"

sum_old_radii= sum_radii;

/*164:*/
#line 237 "./matrixinverse.w"

rad(radx,x);
v_bias::round_nearest();
sum_radii= accumulate(radx.begin(),radx.end(),0.0);


/*:164*/
#line 263 "./matrixinverse.w"



counter++;
}
iterations= counter;





/*:165*/
#line 188 "./matrixinverse.w"


return true;
}


/*:160*/
#line 434 "./matrixinverse.w"

/*166:*/
#line 288 "./matrixinverse.w"

bool MatrixInverse::encloseMatrixInverse(iMatrix&Ainv,const pMatrix&A)
{
bool b= invertMatrix(C,A);
if(b==false)
return false;

/*167:*/
#line 319 "./matrixinverse.w"

assignM(Ci,C);
multMiMp(B,Ci,A);
subFromId(B);


/*:167*//*168:*/
#line 328 "./matrixinverse.w"

v_bias::round_up();
double beta= inf_normM(B);
if(beta>=1)return false;

/*:168*/
#line 295 "./matrixinverse.w"

/*169:*/
#line 338 "./matrixinverse.w"


for(unsigned int i= 0;i<sizeM(A);i++)

{
getColumn(b0,C,i);

bool b= encloseLS(x,A,b0,B,beta);
if(b==false)
{
#ifdef VNODE_DEBUG
printMessage("Could not enclose the solution to a linear system");
#endif
return false;
}
setColumn(Ainv,x,i);
}




/*:169*/
#line 296 "./matrixinverse.w"


#ifdef VNODE_DEBUG
iMatrix B= Ainv;
multMiMp(B,Ainv,A);
int n= sizeM(A);
for(int i= 0;i<n;i++)
for(int j= 0;j<n;j++)
{
interval b= B[i][j];
if(i==j)
assert(v_bias::subseteq(interval(1.0),b));
else
assert(v_bias::subseteq(interval(0.0),b));
}
#endif

return true;
}



/*:166*/
#line 435 "./matrixinverse.w"

/*170:*/
#line 364 "./matrixinverse.w"

bool MatrixInverse::orthogonalInverse(iMatrix&Ainv,const pMatrix&A)
{
transpose(C,A);
/*167:*/
#line 319 "./matrixinverse.w"

assignM(Ci,C);
multMiMp(B,Ci,A);
subFromId(B);


/*:167*//*168:*/
#line 328 "./matrixinverse.w"

v_bias::round_up();
double beta= inf_normM(B);
if(beta>=1)return false;

/*:168*/
#line 368 "./matrixinverse.w"

/*169:*/
#line 338 "./matrixinverse.w"


for(unsigned int i= 0;i<sizeM(A);i++)

{
getColumn(b0,C,i);

bool b= encloseLS(x,A,b0,B,beta);
if(b==false)
{
#ifdef VNODE_DEBUG
printMessage("Could not enclose the solution to a linear system");
#endif
return false;
}
setColumn(Ainv,x,i);
}




/*:169*/
#line 369 "./matrixinverse.w"

return true;
}



/*:170*/
#line 436 "./matrixinverse.w"


}


#line 105 "./vnode.w"




/*:173*/
