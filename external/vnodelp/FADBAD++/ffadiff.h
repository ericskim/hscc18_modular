// Copyright (C) 1996-2006 Ole Stauning (ole.st@uning.dk)
// All rights reserved.

// This code is provided "as is", without any warranty of any kind,
// either expressed or implied, including but not limited to, any implied
// warranty of merchantibility or fitness for any purpose. In no event
// will any party who distributed the code be liable for damages or for
// any claim(s) by any other party, including but not limited to, any
// lost profits, lost monies, lost data or data rendered inaccurate,
// losses sustained by third parties, or any other special, incidental or
// consequential damages arising out of the use or inability to use the
// program, even if the possibility of such damages has been advised
// against. The entire risk as to the quality, the performance, and the
// fitness of the program for any particular purpose lies with the party
// using the code.

// This code, and any derivative of this code, may not be used in a
// commercial package without the prior explicit written permission of
// the authors. Verbatim copies of this code may be made and distributed
// in any medium, provided that this copyright notice is not removed or
// altered in any way. No fees may be charged for distribution of the
// codes, other than a fee to cover the cost of the media and a
// reasonable handling fee.

// ***************************************************************
// ANY USE OF THIS CODE CONSTITUTES ACCEPTANCE OF THE TERMS OF THE
//                         COPYRIGHT NOTICE
// ***************************************************************

#ifndef _FFADIFF_H
#define _FFADIFF_H

#include <iostream>

#include "fadbad.h"

template <class T, int N>
class FFTypeName
{
public:
	T v;
	T g[N];
	bool depend;

	FFTypeName():v(Op<T>::myZero()),depend(false){}
	template <class U> /*explicit*/ FFTypeName(const U& x):v(x),depend(false){}
	FFTypeName(const FFTypeName<T,N>& a);
	FFTypeName<T,N>& operator= (const FFTypeName<T,N>& a);
	template <class U> FFTypeName<T,N>& operator= (const U& a);
	void diff(int idx);
	T& x();
	T& d(int i);

	template <class U> FFTypeName<T,N>& operator += (const U& x);
	template <class U> FFTypeName<T,N>& operator -= (const U& x);
	template <class U> FFTypeName<T,N>& operator *= (const U& x);
	template <class U> FFTypeName<T,N>& operator /= (const U& x);
};

template <class T, int N>
INLINE2 FFTypeName<T,N>::FFTypeName(const FFTypeName<T,N>& a):v(a.v),depend(a.depend)
{
	if (depend) for(int i=0;i<N;i++) g[i]=a.g[i];
}

template <class T, int N>
INLINE2 FFTypeName<T,N>& FFTypeName<T,N>::operator= (const FFTypeName<T,N>& a)
{
	if (this==&a) return *this;
	v=a.v;
	depend=a.depend;
	if (depend) for(int i=0;i<N;i++) g[i]=a.g[i];
	return *this;
}

template <class T, int N> template <class U>
INLINE2 FFTypeName<T,N>& FFTypeName<T,N>::operator= (const U& a)
{
	v=a;
	depend=false;
	return *this;
}

template <class T, int N>
INLINE2 void FFTypeName<T,N>::diff(int idx)
{
	for(int i=0;i<idx;g[i++]=Op<T>::myZero());
	g[idx]=Op<T>::myOne();
	for(int i=idx+1;i<N;g[i++]=Op<T>::myZero());
	depend=true;
}

template <class T, int N>
INLINE1 T& FFTypeName<T,N>::x()
{
	return v;
}

template <class T, int N>
INLINE2 T& FFTypeName<T,N>::d(int i)
{
	USER_ASSERT(i<N,"Index out of bounds");
	if (depend)
	{
		return g[i];
	}
	else
	{
		static T zero(Op<T>::myZero());
		return zero;
	}
}

template <class T, int N>
INLINE2 std::ostream& operator << (std::ostream& os, const FFTypeName<T,N>& o)
{
	os << "{" << o.v;
	if (o.depend)
	{
		os << " [ ";
		for (int i=0;i<N;i++) os << o.g[i] << " ";
		os << "]";
	}
	os << "}";
	return os;
}

template <class T, class U, int N>
INLINE2 FFTypeName<T,N> add1(const U& a, const FFTypeName<T,N>& b)
{
	FFTypeName<T,N> c(a+b.v);
	if (!b.depend) return c;
	for (int i=0;i<N;i++) c.g[i]=b.g[i];
	c.depend=true;
	return c;
}

template <class T, class U, int N>
INLINE2 FFTypeName<T,N> add2(const FFTypeName<T,N>& a,const U& b)
{
	FFTypeName<T,N> c(a.v+b);
	if (!a.depend) return c;
	for (int i=0;i<N;i++) c.g[i]=a.g[i];
	c.depend=true;
	return c;
}

template <class T, int N>
INLINE2 FFTypeName<T,N> operator+ (const typename Op<T>::Base& a, const FFTypeName<T,N>& b)
{
	return add1(a,b);
}
/*
template <class T, int N>
INLINE2 FFTypeName<T,N> operator+ (const T& a, const FFTypeName<T,N>& b)
{
	return add1(a,b);
}
*/
template <class T, int N>
INLINE2 FFTypeName<T,N> operator+ (const FFTypeName<T,N>& a,const typename Op<T>::Base& b)
{
	return add2(a,b);
}
/*
template <class T, int N>
INLINE2 FFTypeName<T,N> operator+ (const FFTypeName<T,N>& a,const T& b)
{
	return add2(a,b);
}
*/
template <class T, int N>
INLINE2 FFTypeName<T,N> operator+ (const FFTypeName<T,N>& a, const FFTypeName<T,N>& b)
{
	if (!a.depend && !b.depend) return FFTypeName<T,N>(a.v+b.v);
	if (!a.depend) return add1(a.v,b);
	if (!b.depend) return add2(a,b.v);
	FFTypeName<T,N> c(a.v+b.v);
	for (int i=0;i<N;i++) c.g[i]=a.g[i]+b.g[i];
	c.depend=true;
	return c;
}

template <class T, class U, int N>
INLINE2 FFTypeName<T,N> sub1(const U& a, const FFTypeName<T,N>& b)
{
	FFTypeName<T,N> c(a-b.v);
	if (!b.depend) return c;
	for (int i=0;i<N;i++) c.g[i]=Op<T>::myNeg(b.g[i]);
	c.depend=true;
	return c;
}

template <class T, class U, int N>
INLINE2 FFTypeName<T,N> sub2(const FFTypeName<T,N>& a,const U& b)
{
	FFTypeName<T,N> c(a.v-b);
	if (!a.depend) return c;
	for (int i=0;i<N;i++) c.g[i]=a.g[i];
	c.depend=true;
	return c;
}

template <class T, int N>
INLINE2 FFTypeName<T,N> operator- (const typename Op<T>::Base& a, const FFTypeName<T,N>& b)
{
	return sub1(a,b);
}
/*
template <class T, int N>
INLINE2 FFTypeName<T,N> operator- (const T& a, const FFTypeName<T,N>& b)
{
	return sub1(a,b);
}
*/
template <class T, int N>
INLINE2 FFTypeName<T,N> operator- (const FFTypeName<T,N>& a,const typename Op<T>::Base& b)
{
	return sub2(a,b);
}
/*
template <class T, int N>
INLINE2 FFTypeName<T,N> operator- (const FFTypeName<T,N>& a,const T& b)
{
	return sub2(a,b);
}
*/
template <class T, int N>
INLINE2 FFTypeName<T,N> operator- (const FFTypeName<T,N>& a, const FFTypeName<T,N>& b)
{
	if (!a.depend && !b.depend) return FFTypeName<T,N>(a.v-b.v);
	if (!a.depend) return sub1(a.v,b);
	if (!b.depend) return sub2(a,b.v);
	FFTypeName<T,N> c(a.v-b.v);
	for (int i=0;i<N;i++) c.g[i]=a.g[i]-b.g[i];
	c.depend=true;
	return c;
}

template <class T, class U, int N>
INLINE2 FFTypeName<T,N> mul1(const U& a, const FFTypeName<T,N>& b)
{
	FFTypeName<T,N> c(a*b.v);
	if (!b.depend) return c;
	for (int i=0;i<N;i++) c.g[i]=b.g[i]*a;
	c.depend=true;
	return c;
}

template <class T, class U, int N>
INLINE2 FFTypeName<T,N> mul2(const FFTypeName<T,N>& a,const U& b)
{
	FFTypeName<T,N> c(a.v*b);
	if (!a.depend) return c;
	for (int i=0;i<N;i++) c.g[i]=a.g[i]*b;
	c.depend=true;
	return c;
}

template <class T, int N>
INLINE2 FFTypeName<T,N> operator* (const typename Op<T>::Base& a, const FFTypeName<T,N>& b)
{
	return mul1(a,b);
}
/*
template <class T, int N>
INLINE2 FFTypeName<T,N> operator* (const T& a, const FFTypeName<T,N>& b)
{
	return mul1(a,b);
}
*/
template <class T, int N>
INLINE2 FFTypeName<T,N> operator* (const FFTypeName<T,N>& a,const typename Op<T>::Base& b)
{
	return mul2(a,b);
}
/*
template <class T, int N>
INLINE2 FFTypeName<T,N> operator* (const FFTypeName<T,N>& a,const T& b)
{
	return mul2(a,b);
}
*/
template <class T, int N>
INLINE2 FFTypeName<T,N> operator* (const FFTypeName<T,N>& a, const FFTypeName<T,N>& b)
{
	if (!a.depend && !b.depend) return FFTypeName<T,N>(a.v*b.v);
	if (!a.depend) return mul1(a.v,b);
	if (!b.depend) return mul2(a,b.v);
	FFTypeName<T,N> c(a.v*b.v);
	for (int i=0;i<N;i++) c.g[i]=a.g[i]*b.v+b.g[i]*a.v;
	c.depend=true;
	return c;
}

template <class T, class U, int N>
INLINE2 FFTypeName<T,N> div1(const U& a, const FFTypeName<T,N>& b)
{
	FFTypeName<T,N> c(a/b.v);
	if (!b.depend) return c;
	T tmp(Op<T>::myNeg(c.v/b.v));
	for (int i=0;i<N;i++) c.g[i]=tmp*b.g[i];
	c.depend=true;
	return c;
}

template <class T, class U, int N>
INLINE2 FFTypeName<T,N> div2(const FFTypeName<T,N>& a, const U& b)
{
	FFTypeName<T,N> c(a.v/b);
	if (!a.depend) return c;
	for (int i=0;i<N;i++) c.g[i]=(a.g[i])/b;
	c.depend=true;
	return c;
}

template <class T, int N>
INLINE2 FFTypeName<T,N> operator/ (const typename Op<T>::Base& a, const FFTypeName<T,N>& b)
{
	return div1(a,b);
}
/*
template <class T, int N>
INLINE2 FFTypeName<T,N> operator/ (const T& a, const FFTypeName<T,N>& b)
{
	return div1(a,b);
}
*/
template <class T, int N>
INLINE2 FFTypeName<T,N> operator/ (const FFTypeName<T,N>& a, const typename Op<T>::Base& b)
{
	return div2(a,b);
}
/*
template <class T, int N>
INLINE2 FFTypeName<T,N> operator/ (const FFTypeName<T,N>& a, const T& b)
{
	return div2(a,b);
}
*/
template <class T, int N>
INLINE2 FFTypeName<T,N> operator/ (const FFTypeName<T,N>& a, const FFTypeName<T,N>& b)
{
	if (!a.depend && !b.depend) return FFTypeName<T,N>(a.v/b.v);
	if (!a.depend) return div1(a.v,b);
	if (!b.depend) return div2(a,b.v);
	FFTypeName<T,N> c(a.v/b.v);
	for (int i=0;i<N;i++) c.g[i]=(a.g[i]-c.v*b.g[i])/b.v;
	c.depend=true;
	return c;
}

template <class T, int N>
INLINE2 FFTypeName<T,N> pow1(const T& a, const FFTypeName<T,N>& b)
{
	FFTypeName<T,N> c(Op<T>::myPow(a,b.v));
	if (!b.depend) return c;
	T tmp(c.v*Op<T>::myLog(a));
	for (int i=0;i<N;i++) c.g[i]=tmp*b.g[i];
	c.depend=true;
	return c;
}
template <class T, int N>
INLINE2 FFTypeName<T,N> pow2(const FFTypeName<T,N>& a, const T& b)
{
	FFTypeName<T,N> c(Op<T>::myPow(a.v,b));
	if (!a.depend) return c;
	T tmp(b*Op<T>::myPow(a.v,b-Op<T>::myOne()));
	for (int i=0;i<N;i++) c.g[i]=tmp*a.g[i];
	c.depend=true;
	return c;
}

template <class T, int N>
INLINE2 FFTypeName<T,N> pow (const typename Op<T>::Base& a, const FFTypeName<T,N>& b)
{
	FFTypeName<T,N> c(Op<T>::myPow(a,b.v));
	if (!b.depend) return c;
	T tmp(c.v*Op<T>::myLog(a));
	for (int i=0;i<N;i++) c.g[i]=tmp*b.g[i];
	c.depend=true;
	return c;
}

template <class T, int N>
INLINE2 FFTypeName<T,N> pow (const FFTypeName<T,N>& a,const typename Op<T>::Base& b)
{
	FFTypeName<T,N> c(Op<T>::myPow(a.v,b));
	if (!a.depend) return c;
	T tmp(b*Op<T>::myPow(a.v,b-Op<T>::myOne()));
	for (int i=0;i<N;i++) c.g[i]=tmp*a.g[i];
	c.depend=true;
	return c;
}

template <class T, int N>
INLINE2 FFTypeName<T,N> pow (const FFTypeName<T,N>& a, const FFTypeName<T,N>& b)
{
	if (!a.depend && !b.depend) return FFTypeName<T,N>(Op<T>::myPow(a.v,b.v));
	if (!a.depend) return pow1(a.v,b);
	if (!b.depend) return pow2(a,b.v);
	FFTypeName<T,N> c(Op<T>::myPow(a.v,b.v));
	T tmp(b.v*Op<T>::myPow(a.v,b.v-Op<T>::myOne())),tmp1(c.v*Op<T>::myLog(a.v));
	for (int i=0;i<N;i++) c.g[i]=tmp*a.g[i]+tmp1*b.g[i];
	c.depend=true;
	return c;
}

/* Unary operators */
template <class T, int N>
INLINE2 FFTypeName<T,N> operator+ (const FFTypeName<T,N>& a)
{
	FFTypeName<T,N> c(a.v);
	if (!a.depend) return c;
	for (int i=0;i<N;i++) c.g[i]=a.g[i];
	c.depend=true;
	return c;
}

template <class T, int N>
INLINE2 FFTypeName<T,N> operator- (const FFTypeName<T,N>& a)
{
	FFTypeName<T,N> c(Op<T>::myNeg(a.v));
	if (!a.depend) return c;
	for (int i=0;i<N;i++) c.g[i]=Op<T>::myNeg(a.g[i]);
	c.depend=true;
	return c;
}

template <class T, int N>
INLINE2 FFTypeName<T,N> pow (const FFTypeName<T,N>& a,int b)
{
	FFTypeName<T,N> c(Op<T>::myPow(a.v,b));
	if (!a.depend) return c;
	T tmp(Op<T>::myInteger(b)*Op<T>::myPow(a.v,b-1));
	for (int i=0;i<N;i++) c.g[i]=a.g[i]*tmp;
	c.depend=true;
	return c;
}

template <class T, int N>
INLINE2 FFTypeName<T,N> sqr (const FFTypeName<T,N>& a)
{
	FFTypeName<T,N> c(Op<T>::mySqr(a.v));
	if (!a.depend) return c;
	T tmp(Op<T>::myTwo()*a.v);
	for (int i=0;i<N;i++) c.g[i]=a.g[i]*tmp;
	c.depend=true;
	return c;
}

template <class T, int N>
INLINE2 FFTypeName<T,N> exp (const FFTypeName<T,N>& a)
{
	FFTypeName<T,N> c(Op<T>::myExp(a.v));
	if (!a.depend) return c;
	for (int i=0;i<N;i++) c.g[i]=a.g[i]*c.v;
	c.depend=true;
	return c;
}

template <class T, int N>
INLINE2 FFTypeName<T,N> log (const FFTypeName<T,N>& a)
{
	FFTypeName<T,N> c(Op<T>::myLog(a.v));
	if (!a.depend) return c;
	for (int i=0;i<N;i++) c.g[i]=a.g[i]/a.v;
	c.depend=true;
	return c;
}

template <class T, int N>
INLINE2 FFTypeName<T,N> sqrt (const FFTypeName<T,N>& a)
{
	FFTypeName<T,N> c(Op<T>::mySqrt(a.v));
	if (!a.depend) return c;
	T tmp(c.v*Op<T>::myTwo());
	for (int i=0;i<N;i++) c.g[i]=a.g[i]/tmp;
	c.depend=true;
	return c;
}

template <class T, int N>
INLINE2 FFTypeName<T,N> sin (const FFTypeName<T,N>& a)
{
	FFTypeName<T,N> c(Op<T>::mySin(a.v));
	if (!a.depend) return c;
	T tmp(Op<T>::myCos(a.v));
	for (int i=0;i<N;i++) c.g[i]=a.g[i]*tmp;
	c.depend=true;
	return c;
}

template <class T, int N>
INLINE2 FFTypeName<T,N> cos (const FFTypeName<T,N>& a)
{
	FFTypeName<T,N> c(Op<T>::myCos(a.v));
	if (!a.depend) return c;
	T tmp(-Op<T>::mySin(a.v));
	for (int i=0;i<N;i++) c.g[i]=a.g[i]*tmp;
	c.depend=true;
	return c;
}

template <class T, int N>
INLINE2 FFTypeName<T,N> tan (const FFTypeName<T,N>& a)
{
	FFTypeName<T,N> c(Op<T>::myTan(a.v));
	if (!a.depend) return c;
	T tmp(Op<T>::myOne()+Op<T>::mySqr(c.v));
	for (int i=0;i<N;i++)  c.g[i]=a.g[i]*tmp;
	c.depend=true;
	return c;
}

template <class T, int N>
INLINE2 FFTypeName<T,N> asin (const FFTypeName<T,N>& a)
{
	FFTypeName<T,N> c(Op<T>::myAsin(a.v));
	if (!a.depend) return c;
	T tmp(Op<T>::myInv(Op<T>::mySqrt(Op<T>::myOne()-Op<T>::mySqr(a.v))));
	for (int i=0;i<N;i++) c.g[i]=a.g[i]*tmp;
	c.depend=true;
	return c;
}

template <class T, int N>
INLINE2 FFTypeName<T,N> acos (const FFTypeName<T,N>& a)
{
	FFTypeName<T,N> c(Op<T>::myAcos(a.v));
	if (!a.depend) return c;
	T tmp(Op<T>::myNeg(Op<T>::myInv(Op<T>::mySqrt(Op<T>::myOne()-Op<T>::mySqr(a.v)))));
	for (int i=0;i<N;i++) c.g[i]=a.g[i]*tmp;
	c.depend=true;
	return c;
}

template <class T, int N>
INLINE2 FFTypeName<T,N> atan (const FFTypeName<T,N>& a)
{
	FFTypeName<T,N> c(Op<T>::myAtan(a.v));
	if (!a.depend) return c;
	T tmp(Op<T>::myInv(Op<T>::myOne()+Op<T>::mySqr(a.v)));
	for (int i=0;i<N;i++) c.g[i]=a.g[i]*tmp;
	c.depend=true;
	return c;
}

/* --------------------------------------------------------------------- */
/*                                                                       */
/*                    COMPOUND ASSIGNMENT OPERATORS                      */
/*                                                                       */
/* --------------------------------------------------------------------- */

template <class T, int N> template <class U>
INLINE1 FFTypeName<T,N>& FFTypeName<T,N>::operator += (const U& x)
{
	return *this = (*this)+x;
}

template <class T, int N> template <class U>
INLINE1 FFTypeName<T,N>& FFTypeName<T,N>::operator -= (const U& x)
{
	return *this = (*this)-x;
}

template <class T, int N> template <class U>
INLINE1 FFTypeName<T,N>& FFTypeName<T,N>::operator *= (const U& x)
{
	return *this = (*this)*x;
}

template <class T, int N> template <class U>
INLINE1 FFTypeName<T,N>& FFTypeName<T,N>::operator /= (const U& x)
{
	return *this = (*this)/x;
}

/* --------------------------------------------------------------------- */
/*                                                                       */
/*                          COMPARISON OPERATORS                         */
/*                                                                       */
/* --------------------------------------------------------------------- */

template <class T, int N>
INLINE1 bool operator == (const FFTypeName<T,N> &a, const FFTypeName<T,N> &b)
{
	return eq(a.v,b.v);
}
template <class T, int N>
INLINE1 bool operator == (const FFTypeName<T,N> &a, const typename Op<T>::Base &b)
{
	return eq(a.v,b);
}
template <class T, int N>
INLINE1 bool operator == (const typename Op<T>::Base &a, const FFTypeName<T,N> &b)
{
	return eq(a,b.v);
}
template <class T, int N>
INLINE1 bool operator != (const FFTypeName<T,N> &a, const FFTypeName<T,N> &b)
{
	return ne(a.v,b.v);
}
template <class T, int N>
INLINE1 bool operator != (const FFTypeName<T,N> &a, const typename Op<T>::Base &b)
{
	return ne(a.v,b);
}
template <class T, int N>
INLINE1 bool operator != (const typename Op<T>::Base &a, const FFTypeName<T,N> &b)
{
	return ne(a,b.v);
}
template <class T, int N>
INLINE1 bool operator > (const FFTypeName<T,N> &a, const FFTypeName<T,N> &b)
{
	return gt(a.v,b.v);
}
template <class T, int N>
INLINE1 bool operator > (const FFTypeName<T,N> &a, const typename Op<T>::Base &b)
{
	return gt(a.v,b);
}
template <class T, int N>
INLINE1 bool operator > (const typename Op<T>::Base &a, const FFTypeName<T,N> &b)
{
	return gt(a,b.v);
}
template <class T, int N>
INLINE1 bool operator >= (const FFTypeName<T,N> &a, const FFTypeName<T,N> &b)
{
	return ge(a.v,b.v);
}
template <class T, int N>
INLINE1 bool operator >= (const FFTypeName<T,N> &a, const typename Op<T>::Base &b)
{
	return ge(a.v,b);
}
template <class T, int N>
INLINE1 bool operator >= (const typename Op<T>::Base &a, const FFTypeName<T,N> &b)
{
	return ge(a,b.v);
}
template <class T, int N>
INLINE1 bool operator < (const FFTypeName<T,N> &a, const FFTypeName<T,N> &b)
{
	return lt(a.v,b.v);
}
template <class T, int N>
INLINE1 bool operator < (const FFTypeName<T,N> &a, const typename Op<T>::Base &b)
{
	return lt(a.v,b);
}
template <class T, int N>
INLINE1 bool operator < (const typename Op<T>::Base &a, const FFTypeName<T,N> &b)
{
	return lt(a,b.v);
}
template <class T, int N>
INLINE1 bool operator <= (const FFTypeName<T,N> &a, const FFTypeName<T,N> &b)
{
	return le(a.v,b.v);
}
template <class T, int N>
INLINE1 bool operator <= (const FFTypeName<T,N> &a, const typename Op<T>::Base &b)
{
	return le(a.v,b);
}
template <class T, int N>
INLINE1 bool operator <= (const typename Op<T>::Base &a, const FFTypeName<T,N> &b)
{
	return le(a,b.v);
}

template <typename U, int N> struct Op< FFTypeName<U,N> >
{
	typedef FFTypeName<U,N> T;
	typedef typename Op<U>::Base Base;
	static Base myInteger(const int i) { return Base(i); }
	static Base myZero() { return myInteger(0); }
	static Base myOne() { return myInteger(1);}
	static Base myTwo() { return myInteger(2); }
	static T myPos(const T& x) { return +x; }
	static T myNeg(const T& x) { return -x; }
	static T& myCadd(T& x, const T& y) { return x+=y; }
	static T& myCsub(T& x, const T& y) { return x-=y; }
	static T& myCmul(T& x, const T& y) { return x*=y; }
	static T& myCdiv(T& x, const T& y) { return x/=y; }
	static T myInv(const T& x) { return myOne()/x; }
	static T mySqr(const T& x) { return x*x; }
	static T myPow(const T& x, const int n) { return ::pow(x,n); }
	static T myPow(const T& x, const T& y) { return ::pow(x,y); }
	static T mySqrt(const T& x) { return ::sqrt(x); }
	static T myLog(const T& x) { return ::log(x); }
	static T myExp(const T& x) { return ::exp(x); }
	static T mySin(const T& x) { return ::sin(x); }
	static T myCos(const T& x) { return ::cos(x); }
	static T myTan(const T& x) { return ::tan(x); }
	static T myAsin(const T& x) { return ::asin(x); }
	static T myAcos(const T& x) { return ::acos(x); }
	static T myAtan(const T& x) { return ::atan(x); }
	static bool myEq(const T& x, const T& y) { return x==y; }
	static bool myNe(const T& x, const T& y) { return x!=y; }
	static bool myLt(const T& x, const T& y) { return x<y; }
	static bool myLe(const T& x, const T& y) { return x<=y; }
	static bool myGt(const T& x, const T& y) { return x>y; }
	static bool myGe(const T& x, const T& y) { return x>=y; }
};

#endif
