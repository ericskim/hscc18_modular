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

//
// NORMAL DISTRIBUTION FUNCTIONS AND AUTOMATIC DIFFERENTIATION OF THESE

#ifndef _NDFAD_H
#define _NDFAD_H

#include "fadbad.h"

#ifdef _FFADIFF_H
// The fast-forward derivative of the normal distribution
template <class T, int N>
INLINE2 FF<T,N> nd(const FF<T,N>& a)
{
	FF<T,N> c(nd(a.v));
	if (!a.depend) return c;
	for (int i=0;i<N;i++) c.g[i]=-a.g[i]*c.v*a.v;
	c.depend=true;
	return c;
}
// The fast-forward derivative of the cumulative normal distribution
template <class T, int N>
INLINE2 FF<T,N> cnd(const FF<T,N>& a) 
{
	FF<T,N> c(cnd(a.v));
	if (!a.depend) return c;
	T tmp(nd(a.v));
	for (int i=0;i<N;i++) c.g[i]=a.g[i]*tmp;
	c.depend=true;
	return c;
}
// The fast-forward derivative of the inverse cumulative normal distribution
template <class T, int N>
INLINE2 FF<T,N> icnd(const FF<T,N>& a) 
{
	FF<T,N> c(icnd(a.v));
	if (!a.depend) return c;
	T tmp(Op<T>::myInv(nd(c.v)));
	for (int i=0;i<N;i++) c.g[i]=a.g[i]*tmp;
	c.depend=true;
	return c;
}
#endif

#ifdef _FADIFF_H
// The forward derivative of the normal distribution
template <class T>
INLINE2 F<T> nd(const F<T>& a)
{
	F<T> c(nd(a.v));
	c.touchg(a.gsize);
	for (int i=0;i<a.gsize;i++) c.g[i]=-a.g[i]*c.v*a.v;
	return c;
}
// The fast-forward derivative of the cumulative normal distribution
template <class T>
INLINE2 F<T> cnd(const F<T>& a) 
{
	F<T> c(cnd(a.v));
	c.touchg(a.gsize);
	T tmp(nd(a.v));
	for (int i=0;i<a.gsize;i++) c.g[i]=a.g[i]*tmp;
	return c;
}
// The forward derivative of the inverse cumulative normal distribution
template <class T>
INLINE2 F<T> icnd(const F<T>& a) 
{
	F<T> c(icnd(a.v));
	c.touchg(a.gsize);
	T tmp(Op<T>::myInv(nd(c.v)));
	for (int i=0;i<a.gsize;i++) c.g[i]=a.g[i]*tmp;
	return c;
}
#endif

#ifdef _BADIFF_H
// The backward derivative of the normal distribution
template <class T>
class CNDBackwardOp: public UNBackwardOp<T>
{
public:
	CNDBackwardOp(const CNDBackwardOp&);
	CNDBackwardOp(const BackwardOp<T>& a);
	virtual ~CNDBackwardOp();
	virtual BackwardOp<T>* copy() const;
	virtual void propagateop();
};

template <class T>
INLINE1 CNDBackwardOp<T>::CNDBackwardOp(const CNDBackwardOp& o):
	UNBackwardOp<T>(o)
{}

template <class T>
INLINE1 CNDBackwardOp<T>::CNDBackwardOp(const BackwardOp<T>& a):
	UNBackwardOp<T>(cnd(a.m_val->value()),a.copy())
{}

template <class T>
INLINE1 CNDBackwardOp<T>::~CNDBackwardOp()
{}

template <class T>
INLINE2 BackwardOp<T>* CNDBackwardOp<T>::copy() const
{
	return newCopy(CNDBackwardOp<T>,*this);
}

template <class T>
INLINE2 void CNDBackwardOp<T>::propagateop()
{
	T tmp(nd((this->m_o1->m_val)->value()));
	this->m_o1->m_val->pav(tmp,*this->m_val);
}

template <class T>
INLINE1 CNDBackwardOp<T> cnd (const BackwardOp<T>& a)
{
	return CNDBackwardOp<T>(a);
}

// The backward derivative of the inverse normal distribution
template <class T>
class ICNDBackwardOp: public UNBackwardOp<T>
{
public:
	ICNDBackwardOp(const ICNDBackwardOp&);
	ICNDBackwardOp(const BackwardOp<T>& a);
	virtual ~ICNDBackwardOp();
	virtual BackwardOp<T>* copy() const;
	virtual void propagateop();
};

template <class T>
INLINE1 ICNDBackwardOp<T>::ICNDBackwardOp(const ICNDBackwardOp& o):
	UNBackwardOp<T>(o)
{}

template <class T>
INLINE1 ICNDBackwardOp<T>::ICNDBackwardOp(const BackwardOp<T>& a):
	UNBackwardOp<T>(icnd(a.m_val->value()),a.copy())
{}

template <class T>
INLINE1 ICNDBackwardOp<T>::~ICNDBackwardOp()
{}

template <class T>
INLINE2 BackwardOp<T>* ICNDBackwardOp<T>::copy() const
{
	return newCopy(ICNDBackwardOp<T>,*this);
}

template <class T>
INLINE2 void ICNDBackwardOp<T>::propagateop()
{
	T tmp(Op<T>::myInv(nd((this->m_val)->value())));
	this->m_o1->m_val->pav(tmp,*this->m_val);
}

template <class T>
INLINE1 ICNDBackwardOp<T> icnd (const BackwardOp<T>& a)
{
	return ICNDBackwardOp<T>(a);
}

#endif

#ifdef _TADIFF_H

template <typename U, int N>
struct CND : public BinTTypeNameHV<U,N>
{
	CND(const U& val, TTypeNameHV<U,N>* pOp, TTypeNameHV<U,N>* pNd):BinTTypeNameHV<U,N>(val,pOp,pNd){}
	unsigned int eval(const unsigned int k)
	{
		unsigned int l=std::max(this->op1Eval(k),this->op2Eval(k));
		if (0==this->length()) { this->val(0)=cnd(this->op1Val(0)); this->length()=1; }
		for(unsigned int i=this->length();i<l;++i)
		{
			this->val(i)=Op<U>::myZero();
			for(unsigned int j=0;j<i;++j) 
			{
				Op<U>::myCadd(this->val(i), Op<U>::myInteger(i-j) * this->op2Val(j) * this->op1Val(i-j));
			}
			Op<U>::myCdiv(this->val(i),Op<U>::myInteger(i));
		}
		return this->length()=l;
	}
};
template <typename U, int N>
TTypeName<U,N> cnd(const TTypeName<U,N>& val)
{
	TTypeName<U,N> tmp(exp(-0.5*val*val)/sqrt(2*M_PI));
	return TTypeName<U,N>(static_cast<TTypeNameHV<U,N>*>(new CND<U,N>(cnd(val.val()), val.getTTypeNameHV(), tmp.getTTypeNameHV())));
}

template <typename U, int N>
struct ICND : public UnTTypeNameHV<U,N>
{
	U m_SQRME[N]; // The coefficients of (-f²/2)
	U m_EXPVAL[N]; // The coefficients of (exp(-f²/2))
	ICND(const U& val, TTypeNameHV<U,N>* pOp):UnTTypeNameHV<U,N>(val,pOp)
	{
		m_SQRME[0]=Op<U>::mySqr(this->val(0));
		m_EXPVAL[0]=Op<U>::myExp(-0.5*m_SQRME[0]);		
	}
	unsigned int eval(const unsigned int k)
	{
		unsigned int l=this->opEval(k);
		if (0==this->length()) 
		{ 
			this->val(0)=Op<U>::myExp(this->opVal(0)); 
			this->length()=1; 
			m_SQRME[0]=Op<U>::mySqr(this->val(0));
			m_EXPVAL[0]=Op<U>::myExp(-0.5*m_SQRME[0]);
		}
		for(unsigned int i=this->length();i<l;++i)
		{
			this->val(i)=Op<U>::myZero();
			if (i>1)
			{
				unsigned int i1=i-1;
				m_SQRME[i1]=Op<U>::myZero();
				unsigned int m=(i1+1)/2;
				for(unsigned int j=0;j<m;++j)
					Op<U>::myCadd(m_SQRME[i1], this->val(i1-j)*this->val(j));
				Op<U>::myCmul(m_SQRME[i1],Op<U>::myTwo());
				if (0==i1%2) Op<U>::myCadd(m_SQRME[i1], Op<U>::mySqr(this->val(m)));
	
				m_EXPVAL[i1]=Op<U>::myZero();
				for(unsigned int j=0;j<i1;++j)
					Op<U>::myCadd(m_EXPVAL[i1], Op<U>::myInteger(i1-j)*m_SQRME[i1-j]*m_EXPVAL[j]);
				Op<U>::myCdiv(m_EXPVAL[i1],Op<U>::myInteger(-2*i1));
			}

			for(unsigned int j=1;j<i;++j) 
				Op<U>::myCsub(this->val(i), Op<U>::myInteger(j) * m_EXPVAL[i-j] * this->val(j));
			Op<U>::myCdiv(this->val(i), Op<U>::myInteger(i) );
			Op<U>::myCadd(this->val(i), this->opVal(i) * sqrt(2*M_PI) );

			Op<U>::myCdiv(this->val(i), m_EXPVAL[0]);
		}
		return this->length()=l;
	}
};
template <typename U, int N>
TTypeName<U,N> icnd(const TTypeName<U,N>& val)
{ 
	return TTypeName<U,N>(static_cast<TTypeNameHV<U,N>*>(new ICND<U,N>(icnd(val.val()), val.getTTypeNameHV())));
}

#endif


#endif
