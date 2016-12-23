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

#ifndef _FADIFF_H
#define _FADIFF_H

#include <iostream>

#include "fadbad.h"

template <typename T>
class FTypeName
{
public:
	T v;
	T *g;
	int gsize;

	FTypeName();
	template <typename U> /*explicit*/ FTypeName(const U& x):v(x),g(0),gsize(0){}
	FTypeName(const FTypeName<T>& a);
	virtual ~FTypeName();
	FTypeName<T>& operator= (const FTypeName<T>& a);
	template <typename U> FTypeName<T>& operator= (const U& a);
	void diff(int idx, int n);
	void touchg(int n);
	T& x();
	T& d(int i);

	template <typename U> FTypeName<T>& operator += (const U& x);
	template <typename U> FTypeName<T>& operator -= (const U& x);
	template <typename U> FTypeName<T>& operator *= (const U& x);
	template <typename U> FTypeName<T>& operator /= (const U& x);
};

template <typename T>
INLINE1 FTypeName<T>::FTypeName():v(Op<T>::myZero()),g(0),gsize(0)
{
}

template <typename T>
INLINE2 FTypeName<T>::FTypeName(const FTypeName& a)
{
	v=a.v;
	gsize=a.gsize;
	if (gsize>0)
	{
		g=newArray(T,gsize);
		for (int i=0;i<gsize;i++) g[i]=a.g[i];
	}
	else
		g=0;
}

template <typename T>
INLINE1 FTypeName<T>::~FTypeName()
{
	if (g!=0) delArray(T,g);
}

template <typename T>
INLINE2 FTypeName<T>& FTypeName<T>::operator= (const FTypeName<T>& a)
{
	if (this==&a) return *this;
	v=a.v;
	if (g!=0) delArray(T,g);
	gsize=a.gsize;
	if (gsize>0)
	{
		g=newArray(T,gsize);
		for (int i=0;i<gsize;i++) g[i]=a.g[i];
	}
	else
		g=0;
	return *this;
}

template <typename T> template <typename U>
INLINE2 FTypeName<T>& FTypeName<T>::operator= (const U& a)
{
	v=a;
	if (g!=0) delArray(T,g);
	gsize=0;
	g=0;
	return *this;
}

template <typename T>
INLINE2 void FTypeName<T>::diff(int idx, int n)
{
	int i;
	g=newArray(T,n);
	gsize=n;
	for(i=0;i<idx;g[i++]=Op<T>::myZero());
	g[idx]=Op<T>::myOne();
	for(i=idx+1;i<n;g[i++]=Op<T>::myZero());
}

template <typename T>
INLINE2 void FTypeName<T>::touchg(int n)
{
// This routine creates the vector for storage of g and initializes
// all elements to zero.
	INTERNAL_ASSERT(g==0,"g vector already allocated");
	if (n>0)
	{
		g=newArray(T,n);
		gsize=n;
	}
}

template <typename T>
INLINE1 T& FTypeName<T>::x()
{
	return v;
}

template <typename T>
INLINE2 T& FTypeName<T>::d(int i)
{
	if (g==0)
	{
		static T zero(Op<T>::myZero());
		return zero;
	}
	else
	{
		USER_ASSERT(i<gsize,"Index out of bounds");
		return g[i];
	}
}

template <typename T>
INLINE2 std::ostream& operator << (std::ostream& os, const FTypeName<T>& o)
{
	os << "{" << o.v;
	if (o.g!=0)
	{
		os << " [ ";
		for (int i=0;i<o.gsize;i++) os << o.g[i] << " ";
		os << "]";
	}
	os << "}";
	return os;
}

template <typename T, typename U>
INLINE2 FTypeName<T> add1(const U& a, const FTypeName<T>& b)
{
	FTypeName<T> c(a+b.v);
	c.touchg(b.gsize);
	for (int i=0;i<b.gsize;i++) c.g[i]=b.g[i];
	return c;
}

template <typename T, typename U>
INLINE2 FTypeName<T> add2(const FTypeName<T>& a,const U& b)
{
	FTypeName<T> c(a.v+b);
	c.touchg(a.gsize);
	for (int i=0;i<a.gsize;i++) c.g[i]=a.g[i];
	return c;
}
/*
template <typename T>
INLINE2 FTypeName<T> operator+ (const typename Op<T>::Base& a, const FTypeName<T>& b)
{
	return add1(a,b);
}
*/

template <typename T, typename U>
INLINE2 FTypeName<T> operator+ (const U& a, const FTypeName<T>& b)
{
	return add1(a,b);
}
/*
template <typename T>
INLINE2 FTypeName<T> operator+ (const FTypeName<T>& a, const typename Op<T>::Base& b)
{
	return add2(a,b);
}
*/

template <typename T, typename U>
INLINE2 FTypeName<T> operator+ (const FTypeName<T>& a, const U& b)
{
	return add2(a,b);
}

template <typename T>
INLINE2 FTypeName<T> operator+ (const FTypeName<T>& a, const FTypeName<T>& b)
{
	if (a.gsize==0) return add1(a.v,b);
	if (b.gsize==0) return add2(a,b.v);
	FTypeName<T> c(a.v+b.v);
	USER_ASSERT(a.gsize==b.gsize,"derivative vectors not of same size in operator+");
	c.touchg(a.gsize);
	for (int i=0;i<a.gsize;i++) c.g[i]=a.g[i]+b.g[i];
	return c;
}

template <typename T, typename U>
INLINE2 FTypeName<T> sub1(const U& a, const FTypeName<T>& b)
{
	FTypeName<T> c(a-b.v);
	c.touchg(b.gsize);
	for (int i=0;i<b.gsize;i++) c.g[i]=Op<T>::myNeg(b.g[i]);
	return c;
}

template <typename T, typename U>
INLINE2 FTypeName<T> sub2(const FTypeName<T>& a,const U& b)
{
	FTypeName<T> c(a.v-b);
	c.touchg(a.gsize);
	for (int i=0;i<a.gsize;i++) c.g[i]=a.g[i];
	return c;
}
/*
template <typename T>
INLINE2 FTypeName<T> operator- (const typename Op<T>::Base& a, const FTypeName<T>& b)
{
	return sub1(a,b);
}
*/

template <typename T, typename U>
INLINE2 FTypeName<T> operator- (const U& a, const FTypeName<T>& b)
{
	return sub1(a,b);
}
/*
template <typename T>
INLINE2 FTypeName<T> operator- (const FTypeName<T>& a,const typename Op<T>::Base& b)
{
	return sub2(a,b);
}
*/

template <typename T, typename U>
INLINE2 FTypeName<T> operator- (const FTypeName<T>& a,const U& b)
{
	return sub2(a,b);
}

template <typename T>
INLINE2 FTypeName<T> operator- (const FTypeName<T>& a, const FTypeName<T>& b)
{
	if (a.gsize==0) return sub1(a.v,b);
	if (b.gsize==0) return sub2(a,b.v);
	FTypeName<T> c(a.v-b.v);
	USER_ASSERT(a.gsize==b.gsize,"derivative vectors not of same size in operator-");
	c.touchg(a.gsize);
	for (int i=0;i<a.gsize;i++) c.g[i]=a.g[i]-b.g[i];
	return c;
}

template <typename T, typename U>
INLINE2 FTypeName<T> mul1(const U& a, const FTypeName<T>& b)
{
	FTypeName<T> c(a*b.v);
	c.touchg(b.gsize);
	for (int i=0;i<b.gsize;i++) c.g[i]=b.g[i]*a;
	return c;
}

template <typename T, typename U>
INLINE2 FTypeName<T> mul2(const FTypeName<T>& a,const U& b)
{
	FTypeName<T> c(a.v*b);
	c.touchg(a.gsize);
	for (int i=0;i<a.gsize;i++) c.g[i]=a.g[i]*b;
	return c;
}
/*
template <typename T>
INLINE2 FTypeName<T> operator* (const typename Op<T>::Base& a, const FTypeName<T>& b)
{
	return mul1(a,b);
}
*/

template <typename T, typename U>
INLINE2 FTypeName<T> operator* (const U& a, const FTypeName<T>& b)
{
	return mul1(a,b);
}
/*
template <typename T>
INLINE2 FTypeName<T> operator* (const FTypeName<T>& a, const typename Op<T>::Base& b)
{
	return mul2(a,b);
}
*/

template <typename T, typename U>
INLINE2 FTypeName<T> operator* (const FTypeName<T>& a, const U& b)
{
	return mul2(a,b);
}

template <typename T>
INLINE2 FTypeName<T> operator* (const FTypeName<T>& a, const FTypeName<T>& b)
{
	if (a.gsize==0) return mul1(a.v,b);
	if (b.gsize==0) return mul2(a,b.v);
	FTypeName<T> c(a.v*b.v);
	USER_ASSERT(a.gsize==b.gsize,"derivative vectors not of same size in operator*");
	c.touchg(a.gsize);
	for (int i=0;i<a.gsize;i++) c.g[i]=a.g[i]*b.v+b.g[i]*a.v;
	return c;
}

template <typename T, typename U>
INLINE2 FTypeName<T> div1(const U& a, const FTypeName<T>& b)
{
	FTypeName<T> c(a/b.v);
	T tmp(Op<T>::myNeg(c.v/b.v));
	c.touchg(b.gsize);
	for (int i=0;i<b.gsize;i++) c.g[i]=tmp*b.g[i];
	return c;
}

template <typename T, typename U>
INLINE2 FTypeName<T> div2(const FTypeName<T>& a, const U& b)
{
	FTypeName<T> c(a.v/b);
	c.touchg(a.gsize);
	for (int i=0;i<a.gsize;i++) c.g[i]=(a.g[i])/b;
	return c;
}
/*
template <typename T>
INLINE2 FTypeName<T> operator/ (const typename Op<T>::Base& a, const FTypeName<T>& b)
{
	return div1(a,b);
}
*/

template <typename T, typename U>
INLINE2 FTypeName<T> operator/ (const U& a, const FTypeName<T>& b)
{
	return div1(a,b);
}
/*
template <typename T>
INLINE2 FTypeName<T> operator/ (const FTypeName<T>& a, const typename Op<T>::Base& b)
{
	return div2(a,b);
}
*/

template <typename T, typename U>
INLINE2 FTypeName<T> operator/ (const FTypeName<T>& a, const U& b)
{
	return div2(a,b);
}

template <typename T>
INLINE2 FTypeName<T> operator/ (const FTypeName<T>& a, const FTypeName<T>& b)
{
	if (a.gsize==0) return div1(a.v,b);
	if (b.gsize==0) return div2(a,b.v);
	FTypeName<T> c(a.v/b.v);
	USER_ASSERT(a.gsize==b.gsize,"derivative vectors not of same size in operator/");
	c.touchg(a.gsize);
	for (int i=0;i<a.gsize;i++) c.g[i]=(a.g[i]-c.v*b.g[i])/b.v;
	return c;
}

template <typename T, typename U>
INLINE2 FTypeName<T> pow1(const U& a, const FTypeName<T>& b)
{
	FTypeName<T> c(Op<T>::myPow(a,b.v));
	T tmp(c.v*Op<T>::myLog(a));
	c.touchg(b.gsize);
	for (int i=0;i<b.gsize;i++) c.g[i]=tmp*b.g[i];
	return c;
}
template <typename T, typename U>
INLINE2 FTypeName<T> pow2(const FTypeName<T>& a, const U& b)
{
	FTypeName<T> c(Op<T>::myPow(a.v,b));
	T tmp(b*Op<T>::myPow(a.v,b-Op<T>::myOne()));
	c.touchg(a.gsize);
	for (int i=0;i<a.gsize;i++) c.g[i]=tmp*a.g[i];
	return c;
}

template <typename T, typename U>
INLINE2 FTypeName<T> pow (const U& a, const FTypeName<T>& b)
{
	return pow1(a,b);
}

template <typename T, typename U>
INLINE2 FTypeName<T> pow (const FTypeName<T>& a,const U& b)
{
	return pow2(a,b);
}

template <typename T>
INLINE2 FTypeName<T> pow (const FTypeName<T>& a, const FTypeName<T>& b)
{
	if (a.gsize==0) return pow1(a.v,b);
	if (b.gsize==0) return pow2(a,b.v);
	FTypeName<T> c(pow(a.v,b.v));
	USER_ASSERT(a.gsize==b.gsize,"derivative vectors not of same size in pow");
	T tmp(b.v*Op<T>::myPow(a.v,b.v-Op<T>::myOne())),tmp1(c.v*Op<T>::myLog(a.v));
	c.touchg(a.gsize);
	for (int i=0;i<a.gsize;i++)
		c.g[i]=tmp*a.g[i]+tmp1*b.g[i];
	return c;
}

/* Unary operators */
template <typename T>
INLINE2 FTypeName<T> operator+ (const FTypeName<T>& a)
{
	FTypeName<T> c(a.v);
	c.touchg(a.gsize);
	for (int i=0;i<a.gsize;i++) c.g[i]=a.g[i];
	return c;
}

template <typename T>
INLINE2 FTypeName<T> operator- (const FTypeName<T>& a)
{
	FTypeName<T> c(Op<T>::myNeg(a.v));
	c.touchg(a.gsize);
	for (int i=0;i<a.gsize;i++) c.g[i]=Op<T>::myNeg(a.g[i]);
	return c;
}

template <typename T>
INLINE2 FTypeName<T> pow (const FTypeName<T>& a,int b)
{
	FTypeName<T> c(Op<T>::myPow(a.v,b));
	c.touchg(a.gsize);
	T tmp(Op<T>::myInteger(b)*Op<T>::myPow(a.v,b-1));
	for (int i=0;i<a.gsize;i++) c.g[i]=a.g[i]*tmp;
	return c;
}

template <typename T>
INLINE2 FTypeName<T> sqr (const FTypeName<T>& a)
{
	FTypeName<T> c(Op<T>::mySqr(a.v));
	c.touchg(a.gsize);
	T tmp(Op<T>::myTwo()*a.v);
	for (int i=0;i<a.gsize;i++) c.g[i]=a.g[i]*tmp;
	return c;
}

template <typename T>
INLINE2 FTypeName<T> exp (const FTypeName<T>& a)
{
	FTypeName<T> c(Op<T>::myExp(a.v));
	c.touchg(a.gsize);
	for (int i=0;i<a.gsize;i++) c.g[i]=a.g[i]*c.v;
	return c;
}

template <typename T>
INLINE2 FTypeName<T> log (const FTypeName<T>& a)
{
	FTypeName<T> c(Op<T>::myLog(a.v));
	c.touchg(a.gsize);
	for (int i=0;i<a.gsize;i++) c.g[i]=a.g[i]/a.v;
	return c;
}

template <typename T>
INLINE2 FTypeName<T> sqrt (const FTypeName<T>& a)
{
	FTypeName<T> c(Op<T>::mySqrt(a.v));
	T tmp(c.v*Op<T>::myTwo());
	c.touchg(a.gsize);
	for (int i=0;i<a.gsize;i++) c.g[i]=a.g[i]/tmp;
	return c;
}

template <typename T>
INLINE2 FTypeName<T> sin (const FTypeName<T>& a)
{
	FTypeName<T> c(Op<T>::mySin(a.v));
	T tmp(Op<T>::myCos(a.v));
	c.touchg(a.gsize);
	for (int i=0;i<a.gsize;i++) c.g[i]=a.g[i]*tmp;
	return c;
}

template <typename T>
INLINE2 FTypeName<T> cos (const FTypeName<T>& a)
{
	FTypeName<T> c(Op<T>::myCos(a.v));
	T tmp(-Op<T>::mySin(a.v));
	c.touchg(a.gsize);
	for (int i=0;i<a.gsize;i++) c.g[i]=a.g[i]*tmp;
	return c;
}

template <typename T>
INLINE2 FTypeName<T> tan (const FTypeName<T>& a)
{
	FTypeName<T> c(Op<T>::myTan(a.v));
	c.touchg(a.gsize);
	T tmp(Op<T>::myOne()+Op<T>::mySqr(c.v));
	for (int i=0;i<a.gsize;i++)  c.g[i]=a.g[i]*tmp;
	return c;
}

template <typename T>
INLINE2 FTypeName<T> asin (const FTypeName<T>& a)
{
	FTypeName<T> c(Op<T>::myAsin(a.v));
	c.touchg(a.gsize);
	T tmp(Op<T>::myInv(sqrt(Op<T>::myOne()-Op<T>::mySqr(a.v))));
	for (int i=0;i<a.gsize;i++) c.g[i]=a.g[i]*tmp;
	return c;
}

template <typename T>
INLINE2 FTypeName<T> acos (const FTypeName<T>& a)
{
	FTypeName<T> c(Op<T>::myAcos(a.v));
	c.touchg(a.gsize);
	T tmp(Op<T>::myNeg(Op<T>::myInv(Op<T>::mySqrt(Op<T>::myOne()-Op<T>::mySqr(a.v)))));
	for (int i=0;i<a.gsize;i++) c.g[i]=a.g[i]*tmp;
	return c;
}

template <typename T>
INLINE2 FTypeName<T> atan (const FTypeName<T>& a)
{
	FTypeName<T> c(Op<T>::myAtan(a.v));
	c.touchg(a.gsize);
	T tmp(Op<T>::myInv(Op<T>::myOne()+Op<T>::mySqr(a.v)));
	for (int i=0;i<a.gsize;i++) c.g[i]=a.g[i]*tmp;
	return c;
}

/* --------------------------------------------------------------------- */
/*                                                                       */
/*                    COMPOUND ASSIGNMENT OPERATORS                      */
/*                                                                       */
/* --------------------------------------------------------------------- */

template <typename T> template <typename U>
INLINE1 FTypeName<T>& FTypeName<T>::operator += (const U& x)
{
	return *this = (*this)+x;
}

template <typename T> template <typename U>
INLINE1 FTypeName<T>& FTypeName<T>::operator -= (const U& x)
{
	return *this = (*this)-x;
}

template <typename T> template <typename U>
INLINE1 FTypeName<T>& FTypeName<T>::operator *= (const U& x)
{
	return *this = (*this)*x;
}

template <typename T> template <typename U>
INLINE1 FTypeName<T>& FTypeName<T>::operator /= (const U& x)
{
	return *this = (*this)/x;
}

/* --------------------------------------------------------------------- */
/*                                                                       */
/*                            COMPARISON OPERATORS                       */
/*                                                                       */
/* --------------------------------------------------------------------- */

template <typename T>
INLINE1 bool operator == (const FTypeName<T> &a, const FTypeName<T> &b)
{
	return Op<T>::myEq(a.v,b.v);
}
template <typename T>
INLINE1 bool operator == (const FTypeName<T> &a, const typename Op<T>::Base &b)
{
	return Op<T>::myEq(a.v,b);
}
template <typename T>
INLINE1 bool operator == (const typename Op<T>::Base &a, const FTypeName<T> &b)
{
	return Op<T>::myEq(a,b.v);
}
template <typename T>
INLINE1 bool operator != (const FTypeName<T> &a, const FTypeName<T> &b)
{
	return Op<T>::myNe(a.v,b.v);
}
template <typename T>
INLINE1 bool operator != (const FTypeName<T> &a, const typename Op<T>::Base &b)
{
	return Op<T>::myNe(a.v,b);
}
template <typename T>
INLINE1 bool operator != (const typename Op<T>::Base &a, const FTypeName<T> &b)
{
	return Op<T>::myNe(a,b.v);
}
template <typename T>
INLINE1 bool operator > (const FTypeName<T> &a, const FTypeName<T> &b)
{
	return Op<T>::myGt(a.v,b.v);
}
template <typename T>
INLINE1 bool operator > (const FTypeName<T> &a, const typename Op<T>::Base &b)
{
	return Op<T>::myGt(a.v,b);
}
template <typename T>
INLINE1 bool operator > (const typename Op<T>::Base &a, const FTypeName<T> &b)
{
	return Op<T>::myGt(a,b.v);
}
template <typename T>
INLINE1 bool operator >= (const FTypeName<T> &a, const FTypeName<T> &b)
{
	return Op<T>::myGe(a.v,b.v);
}
template <typename T>
INLINE1 bool operator >= (const FTypeName<T> &a, const typename Op<T>::Base &b)
{
	return Op<T>::myGe(a.v,b);
}
template <typename T>
INLINE1 bool operator >= (const typename Op<T>::Base &a, const FTypeName<T> &b)
{
	return Op<T>::myGe(a,b.v);
}
template <typename T>
INLINE1 bool operator < (const FTypeName<T> &a, const FTypeName<T> &b)
{
	return Op<T>::myLt(a.v,b.v);
}
template <typename T>
INLINE1 bool operator < (const FTypeName<T> &a, const typename Op<T>::Base &b)
{
	return Op<T>::myLt(a.v,b);
}
template <typename T>
INLINE1 bool operator < (const typename Op<T>::Base &a, const FTypeName<T> &b)
{
	return Op<T>::myLt(a,b.v);
}
template <typename T>
INLINE1 bool operator <= (const FTypeName<T> &a, const FTypeName<T> &b)
{
	return Op<T>::myLe(a.v,b.v);
}
template <typename T>
INLINE1 bool operator <= (const FTypeName<T> &a, const typename Op<T>::Base &b)
{
	return Op<T>::myLe(a.v,b);
}
template <typename T>
INLINE1 bool operator <= (const typename Op<T>::Base &a, const FTypeName<T> &b)
{
	return Op<T>::myLe(a,b.v);
}

template <typename U> struct Op< FTypeName<U> >
{
	typedef FTypeName<U> T;
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
