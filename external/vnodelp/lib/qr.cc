/*154:*/
#line 64 "./qr.w"

#include "basiclinalg.h"
#include "vnoderound.h"

namespace vnodelp{
using namespace v_blas;
/*153:*/
#line 13 "./qr.w"

extern"C"
{

void dgeqrf_(int*m,int*n,double*A,int*lda,double*tau,double*work,
int*lwork,int*info);

void dorgqr_(int*m,int*n,int*k,double*A,int*lda,double*tau,
double*work,int*lwork,int*info);
}


bool computeQR(pMatrix&Q,const pMatrix&A)
{
int n= sizeM(A);
int m= n;
int lda= n;
int k= n;
int info;


int lwork= 10*n;


double*tau= new double[n];
double*work= new double[lwork];
double*M= new double[n*n];

v_bias::round_nearest();
matrix2pointer(M,A);

dgeqrf_(&m,&n,M,&lda,tau,work,&lwork,&info);
if(info==0)
{
dorgqr_(&m,&n,&k,M,&lda,tau,work,&lwork,&info);
if(info==0)
pointer2matrix(Q,M);
}

delete[]M;
delete[]work;
delete[]tau;

if(info==0)
return true;
return false;
}


/*:153*/
#line 70 "./qr.w"

}
#line 1 "./matrixinverse.w"
/*:154*/
