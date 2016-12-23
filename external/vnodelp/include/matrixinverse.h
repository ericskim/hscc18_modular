/*172:*/
#line 406 "./matrixinverse.w"

#ifndef MATRIXINVERSE_H
#define MATRIXINVERSE_H
#include <numeric> 
#include "vnodeinterval.h"
#include "vnoderound.h"
#include "vector_matrix.h"

namespace v_blas
{
using namespace v_bias;
/*156:*/
#line 10 "./matrixinverse.w"

class MatrixInverse
{
public:

MatrixInverse(int n);
bool invertMatrix(pMatrix&Ainv,const pMatrix&A);
bool encloseMatrixInverse(iMatrix&Ainv,const pMatrix&A);
bool orthogonalInverse(iMatrix&Ainv,const pMatrix&A);

~MatrixInverse();

int iterations;
private:
bool encloseLS(iVector&x,const pMatrix&A,const iVector&b0,
const iMatrix&B,double beta);

pVector radx;
iVector b0,x1,x;
iMatrix B,Ci;
pMatrix C;
double*M;

int*ipiv;
double*work;
int lwork;
};

/*:156*/
#line 417 "./matrixinverse.w"

}
#endif

/*:172*/
