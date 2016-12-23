/*137:*/
#line 353 "./linalg.w"

#ifndef BASICLINALG_H
#define BASICLINALG_H
#include <algorithm> 
#include "vector_matrix.h"

namespace v_blas
{

/*129:*/
#line 14 "./linalg.w"


template<class Vector,class scalar> 
inline void setV(Vector&z,scalar a)
{
fill(z.begin(),z.end(),a);
}

template<class Vector1,class Vector2> 
inline void assignV(Vector1&z,const Vector2&x)
{
for(unsigned int i= 0;i<sizeV(z);i++)z[i]= x[i];
}

template<class Vector,class scalar> 
inline void scaleV(Vector&z,scalar a)
{
for(unsigned int i= 0;i<sizeV(z);i++)z[i]*= a;
}


inline void addViVi(iVector&z,const iVector&x)
{
transform(z.begin(),z.end(),x.begin(),z.begin(),plus<v_bias::interval> ());
}

inline void addViVp(iVector&z,const pVector&x)
{
transform(z.begin(),z.end(),x.begin(),z.begin(),plus<v_bias::interval> ());

}


inline void subViVp(iVector&z,const pVector&x)
{
transform(z.begin(),z.end(),x.begin(),z.begin(),minus<v_bias::interval> ());

}

inline void subViVi(iVector&z,const iVector&x)
{
transform(z.begin(),z.end(),x.begin(),z.begin(),minus<v_bias::interval> ());
}


inline void addViVi(iVector&z,const iVector&x,const iVector&y)
{
transform(x.begin(),x.end(),y.begin(),z.begin(),plus<v_bias::interval> ());
}

inline void addViVp(iVector&z,const iVector&x,const pVector&y)
{

transform(x.begin(),x.end(),y.begin(),z.begin(),plus<v_bias::interval> ());

}


inline void subViVp(iVector&z,const iVector&x,const pVector&y)
{

transform(x.begin(),x.end(),y.begin(),z.begin(),minus<v_bias::interval> ());

}



inline double inf_normV(const iVector&z)
{
double s= 0;
for(unsigned int i= 0;i<sizeV(z);i++)
if(v_bias::mag(z[i])> s)
s= v_bias::mag(z[i]);
return s;
}


inline double inf_normV(const pVector&z)
{
double s= 0;
for(unsigned int i= 0;i<sizeV(z);i++)
if(fabs(z[i])> s)
s= fabs(z[i]);
return s;
}


template<class scalar,class Vector1,class Vector2> 
inline void dot_product(scalar&r,const
Vector1&a,const Vector2&b)
{
r= 0.0;
for(unsigned int i= 0;i<a.size();i++)r+= a[i]*b[i];
}


inline double norm2(const pVector&v)
{
double s;
dot_product(s,v,v);
return sqrt(s);
}







/*:129*/
#line 362 "./linalg.w"

/*131:*/
#line 131 "./linalg.w"



inline void multMiVi(iVector&z,const iMatrix&A,
const iVector&x)
{
for(unsigned int i= 0;i<A.size();i++)dot_product(z[i],A[i],x);
}


inline void multMpVi(iVector&z,const pMatrix&A,const iVector&x)
{
for(unsigned int i= 0;i<A.size();i++)dot_product(z[i],A[i],x);
}

inline void multMiVp(iVector&z,const iMatrix&A,const pVector&x)
{
for(unsigned int i= 0;i<A.size();i++)dot_product(z[i],A[i],x);
}


/*:131*/
#line 363 "./linalg.w"

/*132:*/
#line 156 "./linalg.w"


template<class Matrix,class scalar> 
inline void setM(Matrix&C,scalar a)
{
for(unsigned int i= 0;i<C.size();i++)setV(C[i],a);
}

template<class Matrix> 
inline void setId(Matrix&C)
{
setM(C,0.0);
for(unsigned int i= 0;i<sizeM(C);i++)C[i][i]= 1.0;
}

template<class Matrix1,class Matrix2> 
inline void assignM(Matrix1&C,const Matrix2&A)
{
for(unsigned int i= 0;i<C.size();i++)assignV(C[i],A[i]);
}


template<class Matrix,class scalar> 
inline void scaleM(Matrix&C,scalar a)
{
for(unsigned int i= 0;i<C.size();i++)scaleV(C[i],a);
}

template<class Matrix,class Vector> 
void setColumn(Matrix&C,const Vector&z,unsigned int j)
{
for(unsigned int i= 0;i<sizeM(C);i++)
C[i][j]= z[i];
}


template<class Matrix> 
inline void transpose(Matrix&C,const Matrix&A)
{
for(unsigned int i= 0;i<C.size();i++)setColumn(C,A[i],i);
}

template<class Matrix> 
inline void addId(Matrix&C){
for(unsigned int i= 0;i<sizeM(C);i++)
C[i][i]+= 1.0;
}

template<class Matrix> 
inline void subFromId(Matrix&C){
unsigned int n= sizeM(C);

for(unsigned int i= 0;i<n;i++)
for(unsigned int j= 0;j<n;j++)
C[i][j]= -C[i][j];

for(unsigned int i= 0;i<n;i++)
C[i][i]+= 1.0;
}


inline void addMiMi(iMatrix&C,const iMatrix&A)
{
for(unsigned int i= 0;i<C.size();i++)addViVi(C[i],A[i]);
}



inline void subMiMp(iMatrix&C,const pMatrix&A)
{
for(unsigned int i= 0;i<C.size();i++)subViVp(C[i],A[i]);
}


inline void multMiMi(iMatrix&C,const iMatrix&A,const iMatrix&B)
{
unsigned int n= sizeM(A);
for(unsigned int i= 0;i<n;i++)
for(unsigned int j= 0;j<n;j++)
{
C[i][j]= 0.0;
for(unsigned int k= 0;k<n;k++)C[i][j]+= A[i][k]*B[k][j];
}
}

inline void multMiMp(iMatrix&C,const iMatrix&A,const pMatrix&B)
{
unsigned int n= sizeM(A);
for(unsigned int i= 0;i<n;i++)
for(unsigned int j= 0;j<n;j++)
{
C[i][j]= 0.0;
for(unsigned int k= 0;k<n;k++)C[i][j]+= A[i][k]*B[k][j];
}
}


template<class Matrix> inline double inf_normM(const Matrix&C)
{
unsigned int n= sizeM(C);
double m= 0;
for(unsigned int i= 0;i<n;i++)
{
double s= 0;
for(unsigned int j= 0;j<n;j++)
s+= v_bias::mag(C[i][j]);

if(s> m)m= s;
}
return m;
}




/*:132*/
#line 364 "./linalg.w"

/*134:*/
#line 281 "./linalg.w"



template<class Vector,class Matrix> 
void getColumn(Vector&z,const Matrix&C,unsigned int j)
{
for(unsigned int i= 0;i<sizeM(C);i++)
z[i]= C[i][j];
}





/*:134*/
#line 365 "./linalg.w"

/*135:*/
#line 322 "./linalg.w"

template<class scalar,class Matrix> 
inline void matrix2pointer(scalar*M,const Matrix&C)
{
unsigned int n= sizeM(C);

for(unsigned int j= 0;j<n;j++)
for(unsigned int i= 0;i<n;i++)
M[j*n+i]= C[i][j];
}


/*:135*/
#line 366 "./linalg.w"

/*136:*/
#line 335 "./linalg.w"

template<class Matrix,class scalar> 
inline void pointer2matrix(Matrix&C,const scalar*M)
{
unsigned int n= sizeM(C);

for(unsigned int j= 0;j<n;j++)
for(unsigned int i= 0;i<n;i++)
C[i][j]= M[j*n+i];
}





/*:136*/
#line 367 "./linalg.w"

/*366:*/
#line 5 "./misc.w"

#include <iostream> 
using namespace std;
template<class T> 
void printVector(const T&v,const char*s= 0)
{
if(s)
cout<<s<<" = "<<endl;
for(unsigned int i= 0;i<v_blas::sizeV(v);i++)
cout<<v[i]<<endl;
cout<<endl;
}



/*:366*/
#line 368 "./linalg.w"

}


#endif





/*:137*/
