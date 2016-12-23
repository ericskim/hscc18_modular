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

#ifndef _BADIFF_H
#define _BADIFF_H

#include "fadbad.h"

template <class T>
class BackwardNodeNormal
{
	T m_v;
	T *m_g;
	int m_gsize;
	int m_rc;
public:
	BackwardNodeNormal():m_g(0),m_rc(1),m_gsize(0){}
	template <class U> /*explicit*/ BackwardNodeNormal(const U& x):m_v(x),m_g(0),m_rc(1),m_gsize(0){}
	~BackwardNodeNormal();
	void touchg(int n);
	T& value();
	const T& value() const;
	T& deriv(const int idx);
	const T& deriv(const int idx) const;
	int size() const;
	int rc() const;
	int inc_rc();
	int dec_rc();
	void pv(const BackwardNodeNormal& v);
	void mv(const BackwardNodeNormal& v);
	void pav(const T& a, const BackwardNodeNormal& v);
	void mav(const T& a, const BackwardNodeNormal& v);
};


template <class T>
class BackwardNodeSparse
{
	T m_v;
	T** m_g;
	int m_gsize;
	int m_rc;
public:
	BackwardNodeSparse():m_g(0),m_rc(1),m_gsize(0){}
	template <class U> explicit BackwardNodeSparse(const U& x):m_v(x),m_g(0),m_rc(1),m_gsize(0){}
	~BackwardNodeSparse();
	void touchg(int n);
	T& value();
	T& deriv(const int idx);
	const T& deriv(const int idx) const;
	int size() const;
	int rc() const;
	int inc_rc();
	int dec_rc();
	void pv(const BackwardNodeSparse& v);
	void mv(const BackwardNodeSparse& v);
	void pav(const T& a, const BackwardNodeSparse& v);
	void mav(const T& a, const BackwardNodeSparse& v);
};

template <class T>
class BADConfig
{
public:
#ifdef BADSPARSE
	typedef BackwardNodeSparse<T> BackwardNode;
#else
	typedef BackwardNodeNormal<T> BackwardNode;
#endif
};

template <class T>
class BackwardOp
{
// This class is the baseclass of all operation met during the forward
// sweep. It contains a pointer to a BackwardNode where the result of the operation
// is stored, val.
//

public:
	typedef typename BADConfig<T>::BackwardNode BackwardNode;
	BackwardNode *m_val;

protected:
	BackwardOp();
	BackwardOp(const T &x);
	explicit BackwardOp(const int x);
	BackwardOp(const BackwardOp<T> &x);
public:
	virtual ~BackwardOp();
	void inc_rc() const;
	virtual void dec_rc();
	/* Dummy virtual functions */
	virtual BackwardOp<T>* copy() const;
	virtual void propagate();
	virtual void propagateop();
};

template <class T>
class BINBackwardOp: public BackwardOp<T>
{
// This class represents binary operators. It has two branches, one for each
// of the expresions (operators) it operates on. The pointer o1 refers to
// the left operand and o2 to the right operand.
//

public:
	BackwardOp<T> *m_o1,*m_o2;

protected:
	BINBackwardOp();
	BINBackwardOp(const T& x);
	BINBackwardOp(const T& x, BackwardOp<T>* a, BackwardOp<T>* b);
	BINBackwardOp(const BINBackwardOp<T>& x);
public:
	virtual ~BINBackwardOp();
	virtual void propagate();
};

template <class T>
class UNBackwardOp: public BackwardOp<T>
{
// This class represents unary operators. It has one branch for the
// expresion (operator) it operates on. The pointer o1 refers to this.
//

public:
	BackwardOp<T> *m_o1;

protected:
	UNBackwardOp();
	UNBackwardOp(const T& x);
	explicit UNBackwardOp(const int x);
	UNBackwardOp(const T& x, BackwardOp<T>* a);
	UNBackwardOp(const UNBackwardOp<T>& x);
public:
	virtual ~UNBackwardOp();
	virtual void propagate();
};

template <class T>
class BTypeName : public UNBackwardOp<T>
{
	template <class U, class V>
	struct Assign // The assignment operators allowed when T is not the Base-Class
	{
		Assign(BTypeName<U>& dest, const U& x);
		Assign(BTypeName<U>& dest, const BTypeName<U>& x);
		Assign(BTypeName<U>& dest, const BackwardOp<U>& x);
		Assign(BTypeName<U>& dest, const typename Op<U>::Base& x);
	};
	template <class U>
	struct Assign<U, U> // The assignment operators allowed when T is the Base-Class
	{ 
		Assign(BTypeName<U>& dest, const U& x);
		Assign(BTypeName<U>& dest, const BTypeName<U>& x);
		Assign(BTypeName<U>& dest, const BackwardOp<U>& x);
	};

public:
	BTypeName();
	BTypeName(const BTypeName<T>& x);
	BTypeName(const BackwardOp<T>& x);
	explicit BTypeName(const int x);
	template <class U> explicit BTypeName(const U& x);

	virtual ~BTypeName();
	virtual void dec_rc();
	virtual void propagate();
	void diff(int idx, int n);
	virtual BackwardOp<T>* copy() const;
	T& x();
	T& d(int i);

	/* This is all overloadings of the = operator */
	template <class U>
	BTypeName<T>& operator = (const U& x) { Assign< T, typename Op<T>::Base >(*this,x); return *this; }
	// Keep one non-templated assignment operator to prevent the compiler from making the default ones:
	BTypeName<T>& operator = (const BTypeName<T>& x) 
	{ 
		if (this==&x) return *this;
		Assign< T, typename Op<T>::Base >(*this,x); 
		return *this; 
	}

	/* Declaration of overloading of (?)= */
	template <class U> BTypeName<T>& operator += (const U& x);
	template <class U> BTypeName<T>& operator -= (const U& x);
	template <class U> BTypeName<T>& operator *= (const U& x);
	template <class U> BTypeName<T>& operator /= (const U& x);
};


/* --------------------------------------------------------------------- */
/*                                                                       */
/*              BACKWARD DIFFERENTIATION ADDITION FORMULAS               */
/*                                                                       */
/* --------------------------------------------------------------------- */

/* --------------------------- op + op --------------------------------- */

template <class T>
class ADDBackwardOp: public BINBackwardOp<T>
{
public:
	ADDBackwardOp(const ADDBackwardOp<T>&);
	ADDBackwardOp(const BackwardOp<T>& a, const BackwardOp<T>& b);
	virtual ~ADDBackwardOp();
	virtual BackwardOp<T>* copy() const;
	virtual void propagateop();
};

 /* --------------------------- const + op ------------------------------ */

template <class T, class U>
class ADDBackwardOp3: public UNBackwardOp<T>
{
	U m_ax;
public:
	ADDBackwardOp3(const ADDBackwardOp3<T,U>&);
	ADDBackwardOp3(const U& a, const BackwardOp<T>& b);
	virtual ~ADDBackwardOp3();
	virtual BackwardOp<T>* copy() const;
	virtual void propagateop();
};

/* --------------------------- op + const ------------------------------ */

template <class T, class U>
class ADDBackwardOp4: public UNBackwardOp<T>
{
	U m_bx;
public:
	ADDBackwardOp4(const ADDBackwardOp4<T,U>&);
	ADDBackwardOp4(const BackwardOp<T>& a, const U& b);
	virtual ~ADDBackwardOp4();
	virtual BackwardOp<T>* copy() const;
	virtual void propagateop();
};

/* --------------------------------------------------------------------- */
/*                                                                       */
/*            BACKWARD DIFFERENTIATION SUBTRACTION FORMULAS              */
/*                                                                       */
/* --------------------------------------------------------------------- */

/* --------------------------- op - op --------------------------------- */

template <class T>
class SUBBackwardOp: public BINBackwardOp<T>
{
public:
	SUBBackwardOp(const SUBBackwardOp<T>&);
	SUBBackwardOp(const BackwardOp<T>& a, const BackwardOp<T>& b);
	virtual ~SUBBackwardOp();
	virtual BackwardOp<T>* copy() const;
	virtual void propagateop();
};

/* --------------------------- const - op ------------------------------ */

template <class T, class U>
class SUBBackwardOp3: public UNBackwardOp<T>
{
	U m_ax;
public:
	SUBBackwardOp3(const SUBBackwardOp3<T,U>&);
	SUBBackwardOp3(const U& a, const BackwardOp<T>& b);
	virtual ~SUBBackwardOp3();
	virtual BackwardOp<T>* copy() const;
	virtual void propagateop();
};

/* --------------------------- op - const ------------------------------ */

template <class T, class U>
class SUBBackwardOp4: public UNBackwardOp<T>
{
	U m_bx;
public:
	SUBBackwardOp4(const SUBBackwardOp4<T,U>&);
	SUBBackwardOp4(const BackwardOp<T>& a, const U& b);
	virtual ~SUBBackwardOp4();
	virtual BackwardOp<T>* copy() const;
	virtual void propagateop();
};

/* --------------------------------------------------------------------- */
/*                                                                       */
/*          BACKWARD DIFFERENTIATION MULTIPLICATION FORMULAS             */
/*                                                                       */
/* --------------------------------------------------------------------- */

/* --------------------------- op * op --------------------------------- */

template <class T>
class MULBackwardOp: public BINBackwardOp<T>
{
public:
	MULBackwardOp(const MULBackwardOp<T>&);
	MULBackwardOp(const BackwardOp<T>& a, const BackwardOp<T>& b);
	virtual ~MULBackwardOp();
	virtual BackwardOp<T>* copy() const;
	virtual void propagateop();
};

/* --------------------------- const * op ------------------------------ */

template <class T, class U>
class MULBackwardOp3: public UNBackwardOp<T>
{
	U m_ax;
public:
	MULBackwardOp3(const MULBackwardOp3<T,U>&);
	MULBackwardOp3(const U& a, const BackwardOp<T>& b);
	virtual ~MULBackwardOp3();
	virtual BackwardOp<T>* copy() const;
	virtual void propagateop();
};

/* --------------------------- op * const ------------------------------ */

template <class T, class U>
class MULBackwardOp4: public UNBackwardOp<T>
{
	U m_bx;
public:
	MULBackwardOp4(const MULBackwardOp4<T,U>&);
	MULBackwardOp4(const BackwardOp<T>& a, const U& b);
	virtual ~MULBackwardOp4();
	virtual BackwardOp<T>* copy() const;
	virtual void propagateop();
};

/* --------------------------------------------------------------------- */
/*                                                                       */
/*              BACKWARD DIFFERENTIATION DIVISION FORMULAS               */
/*                                                                       */
/* --------------------------------------------------------------------- */

/* --------------------------- op / op --------------------------------- */

template <class T>
class DIVBackwardOp: public BINBackwardOp<T>
{
public:
	DIVBackwardOp(const DIVBackwardOp<T>&);
	DIVBackwardOp(const BackwardOp<T>& a, const BackwardOp<T>& b);
	virtual ~DIVBackwardOp();
	virtual BackwardOp<T>* copy() const;
	virtual void propagateop();
};

/* --------------------------- typename Op<T>::Base / op ------------------------------ */

template <class T, class U>
class DIVBackwardOp3: public UNBackwardOp<T>
{
	U m_ax;
public:
	DIVBackwardOp3(const DIVBackwardOp3<T,U>&);
	DIVBackwardOp3(const U& a, const BackwardOp<T>& b);
	virtual ~DIVBackwardOp3();
	virtual BackwardOp<T>* copy() const;
	virtual void propagateop();
};

/* --------------------------- op / typename Op<T>::Base ------------------------------ */

template <class T, class U>
class DIVBackwardOp4: public UNBackwardOp<T>
{
	U m_bx;
public:
	DIVBackwardOp4(const DIVBackwardOp4<T,U>&);
	DIVBackwardOp4(const BackwardOp<T>& a, const U& b);
	virtual ~DIVBackwardOp4();
	virtual BackwardOp<T>* copy() const;
	virtual void propagateop();
};

/* --------------------------------------------------------------------- */
/*                                                                       */
/*                BACKWARD DIFFERENTIATION POWER FORMULAS                */
/*                                                                       */
/* --------------------------------------------------------------------- */

/* -------------------------- pow(op,BackwardOp) ------------------------------- */

template <class T>
class POWBackwardOp: public BINBackwardOp<T>
{
public:
	POWBackwardOp(const POWBackwardOp&);
	POWBackwardOp(const BackwardOp<T>& a, const BackwardOp<T>& b);
	virtual ~POWBackwardOp();
	virtual BackwardOp<T>* copy() const;
	virtual void propagateop();
};

/* ----------------------- pow(typename Op<T>::Base,BackwardOp) ---------------------------- */

template <class T>
class POWBackwardOp3: public UNBackwardOp<T>
{
	typename Op<T>::Base m_ax;
public:
	POWBackwardOp3(const POWBackwardOp3&);
	POWBackwardOp3(const typename Op<T>::Base& a, const BackwardOp<T>& b);
	virtual ~POWBackwardOp3();
	virtual BackwardOp<T>* copy() const;
	virtual void propagateop();
};

/* ------------------------ pow(op,typename Op<T>::Base) --------------------------- */

template <class T>
class POWBackwardOp4: public UNBackwardOp<T>
{
	typename Op<T>::Base m_bx;
public:
	POWBackwardOp4(const POWBackwardOp4&);
	POWBackwardOp4(const BackwardOp<T>& a, const typename Op<T>::Base& b);
	virtual ~POWBackwardOp4();
	virtual BackwardOp<T>* copy() const;
	virtual void propagateop();
};

/* --------------------------------------------------------------------- */
/*                                                                       */
/*            BACKWARD DIFFERENTIATION UNARY + FORMULA                   */
/*                                                                       */
/* --------------------------------------------------------------------- */

template <class T>
class UADDBackwardOp: public UNBackwardOp<T>
{
public:
	UADDBackwardOp(const UADDBackwardOp&);
	UADDBackwardOp(const BackwardOp<T>& a);
	virtual ~UADDBackwardOp();
	virtual BackwardOp<T>* copy() const;
	virtual void propagateop();
};

/* --------------------------------------------------------------------- */
/*                                                                       */
/*            BACKWARD DIFFERENTIATION UNARY - FORMULA                   */
/*                                                                       */
/* --------------------------------------------------------------------- */

template <class T>
class USUBBackwardOp: public UNBackwardOp<T>
{
public:
	USUBBackwardOp(const USUBBackwardOp&);
	USUBBackwardOp(const BackwardOp<T>& a);
	virtual ~USUBBackwardOp();
	virtual BackwardOp<T>* copy() const;
	virtual void propagateop();
};

/* --------------------------------------------------------------------- */
/*                                                                       */
/*          BACKWARD DIFFERENTIATION pow(op,int) FORMULA                 */
/*                                                                       */
/* --------------------------------------------------------------------- */

template <class T>
class POWopN: public UNBackwardOp<T>
{
	int m_b;
public:
	POWopN(const POWopN&);
	POWopN(const BackwardOp<T>& a, int ib);
	virtual ~POWopN();
	virtual BackwardOp<T>* copy() const;
	virtual void propagateop();
};

/* --------------------------------------------------------------------- */
/*                                                                       */
/*              BACKWARD DIFFERENTIATION sqr FORMULA                     */
/*                                                                       */
/* --------------------------------------------------------------------- */

template <class T>
class SQRBackwardOp: public UNBackwardOp<T>
{
public:
	SQRBackwardOp(const SQRBackwardOp&);
	SQRBackwardOp(const BackwardOp<T>& a);
	virtual ~SQRBackwardOp();
	virtual BackwardOp<T>* copy() const;
	virtual void propagateop();
};

/* --------------------------------------------------------------------- */
/*                                                                       */
/*              BACKWARD DIFFERENTIATION exp FORMULA                     */
/*                                                                       */
/* --------------------------------------------------------------------- */

template <class T>
class EXPBackwardOp: public UNBackwardOp<T>
{
public:
	EXPBackwardOp(const EXPBackwardOp&);
	EXPBackwardOp(const BackwardOp<T>& a);
	virtual ~EXPBackwardOp();
	virtual BackwardOp<T>* copy() const;
	virtual void propagateop();
};

/* --------------------------------------------------------------------- */
/*                                                                       */
/*              BACKWARD DIFFERENTIATION log FORMULA                     */
/*                                                                       */
/* --------------------------------------------------------------------- */

template <class T>
class LOGBackwardOp: public UNBackwardOp<T>
{
public:
	LOGBackwardOp(const LOGBackwardOp&);
	LOGBackwardOp(const BackwardOp<T>& a);
	virtual ~LOGBackwardOp();
	virtual BackwardOp<T>* copy() const;
	virtual void propagateop();
};

/* --------------------------------------------------------------------- */
/*                                                                       */
/*              BACKWARD DIFFERENTIATION sqrt FORMULA                    */
/*                                                                       */
/* --------------------------------------------------------------------- */

template <class T>
class SQRTBackwardOp: public UNBackwardOp<T>
{
public:
	SQRTBackwardOp(const SQRTBackwardOp&);
	SQRTBackwardOp(const BackwardOp<T>& a);
	virtual ~SQRTBackwardOp();
	virtual BackwardOp<T>* copy() const;
	virtual void propagateop();
};

/* --------------------------------------------------------------------- */
/*                                                                       */
/*              BACKWARD DIFFERENTIATION sin FORMULA                     */
/*                                                                       */
/* --------------------------------------------------------------------- */

template <class T>
class SINBackwardOp: public UNBackwardOp<T>
{
public:
	SINBackwardOp(const SINBackwardOp&);
	SINBackwardOp(const BackwardOp<T>& a);
	virtual ~SINBackwardOp();
	virtual BackwardOp<T>* copy() const;
	virtual void propagateop();
};

/* --------------------------------------------------------------------- */
/*                                                                       */
/*              BACKWARD DIFFERENTIATION cos FORMULA                     */
/*                                                                       */
/* --------------------------------------------------------------------- */

template <class T>
class COSBackwardOp: public UNBackwardOp<T>
{
public:
	COSBackwardOp(const COSBackwardOp&);
	COSBackwardOp(const BackwardOp<T>& a);
	virtual ~COSBackwardOp();
	virtual BackwardOp<T>* copy() const;
	virtual void propagateop();
};

/* --------------------------------------------------------------------- */
/*                                                                       */
/*              BACKWARD DIFFERENTIATION tan FORMULA                     */
/*                                                                       */
/* --------------------------------------------------------------------- */

template <class T>
class TANBackwardOp: public UNBackwardOp<T>
{
public:
	TANBackwardOp(const TANBackwardOp&);
	TANBackwardOp(const BackwardOp<T>& a);
	virtual ~TANBackwardOp();
	virtual BackwardOp<T>* copy() const;
	virtual void propagateop();
};

/* --------------------------------------------------------------------- */
/*                                                                       */
/*              BACKWARD DIFFERENTIATION asin FORMULA                    */
/*                                                                       */
/* --------------------------------------------------------------------- */

template <class T>
class ASINBackwardOp: public UNBackwardOp<T>
{
public:
	ASINBackwardOp(const ASINBackwardOp&);
	ASINBackwardOp(const BackwardOp<T>& a);
	virtual ~ASINBackwardOp();
	virtual BackwardOp<T>* copy() const;
	virtual void propagateop();
};

/* --------------------------------------------------------------------- */
/*                                                                       */
/*              BACKWARD DIFFERENTIATION acos FORMULA                    */
/*                                                                       */
/* --------------------------------------------------------------------- */

template <class T>
class ACOSBackwardOp: public UNBackwardOp<T>
{
public:
	ACOSBackwardOp(const ACOSBackwardOp&);
	ACOSBackwardOp(const BackwardOp<T>& a);
	virtual ~ACOSBackwardOp();
	virtual BackwardOp<T>* copy() const;
	virtual void propagateop();
};

/* --------------------------------------------------------------------- */
/*                                                                       */
/*              BACKWARD DIFFERENTIATION atan FORMULA                    */
/*                                                                       */
/* --------------------------------------------------------------------- */

template <class T>
class ATANBackwardOp: public UNBackwardOp<T>
{
public:
	ATANBackwardOp(const ATANBackwardOp&);
	ATANBackwardOp(const BackwardOp<T>& a);
	virtual ~ATANBackwardOp();
	virtual BackwardOp<T>* copy() const;
	virtual void propagateop();
};


// NON-SPARSE BACKWARD NODE (THE DEFAULT)

template <class T>
INLINE2 BackwardNodeNormal<T>::~BackwardNodeNormal()
{
	ASSERT(m_gsize?m_g!=0:m_g==0)
	if (m_g) delArray(T,m_g);
	m_g=0;
}
template <class T>
INLINE2 void BackwardNodeNormal<T>::touchg(int n)
{
	if (!m_g && n>0)
	{
		m_g=newArray(T,n);
		m_gsize=n;
		for(int i=0;i<n;m_g[i++]=Op<T>::myZero());
	}
}
template <class T>
INLINE1 T& BackwardNodeNormal<T>::value(){ return m_v; }
template <class T>
INLINE1 const T& BackwardNodeNormal<T>::value() const { return m_v; }
template <class T>
INLINE1 T& BackwardNodeNormal<T>::deriv(const int idx){ return m_g[idx]; }
template <class T>
INLINE1 const T& BackwardNodeNormal<T>::deriv(const int idx) const { return m_g[idx]; }
template <class T>
INLINE1 int BackwardNodeNormal<T>::size() const { return m_gsize; }
template <class T>
INLINE1 int BackwardNodeNormal<T>::rc() const { return m_rc; }
template <class T>
INLINE1 int BackwardNodeNormal<T>::inc_rc() { return ++m_rc; }
template <class T>
INLINE1 int BackwardNodeNormal<T>::dec_rc() { return --m_rc; }
template <class T>
INLINE2 void BackwardNodeNormal<T>::pv(const BackwardNodeNormal& v)
{
	if (0==v.m_g) return;
	if (0==m_g)
	{
		m_gsize=v.size();
		m_g=newArray(T,m_gsize);
		for(int i=0;i<m_gsize;++i) m_g[i]=v.m_g[i];
	}
	else
	{
		USER_ASSERT(m_gsize==v.size(), ".diff(i,n) called with two different values of"
		" n in same expression; "<< m_gsize << " and " << v.size())
		for(int i=0;i<m_gsize;++i) Op<T>::myCadd(m_g[i],v.m_g[i]);
	}
}
template <class T>
INLINE2 void BackwardNodeNormal<T>::mv(const BackwardNodeNormal& v)
{
	if (0==v.m_g) return;
	if (0==m_g)
	{
		m_gsize=v.size();
		m_g=newArray(T,m_gsize);
		for(int i=0;i<m_gsize;++i) m_g[i]=Op<T>::myNeg(v.m_g[i]);
	}
	else
	{
		USER_ASSERT(m_gsize==v.size(), ".diff(i,n) called with two different values of"
		" n in same expression; "<< m_gsize << " and " << v.size())
		for(int i=0;i<m_gsize;++i) Op<T>::myCsub(m_g[i],v.m_g[i]);
	}
}
template <class T>
INLINE2 void BackwardNodeNormal<T>::pav(const T& a, const BackwardNodeNormal& v)
{
	if (0==v.m_g) return;
	if (0==m_g)
	{
		m_gsize=v.size();
		m_g=newArray(T,m_gsize);
		for(int i=0;i<m_gsize;++i) m_g[i]=a*v.m_g[i];
	}
	else
	{
		USER_ASSERT(m_gsize==v.size(), ".diff(i,n) called with two different values of"
		" n in same expression; "<< m_gsize << " and " << v.size())
		for(int i=0;i<m_gsize;++i) Op<T>::myCadd(m_g[i],a*v.m_g[i]);
	}
}
template <class T>
INLINE2 void BackwardNodeNormal<T>::mav(const T& a, const BackwardNodeNormal& v)
{
	if (0==v.m_g) return;
	if (0==m_g)
	{
		m_gsize=v.size();
		m_g=newArray(T,m_gsize);
		for(int i=0;i<m_gsize;++i) m_g[i]=Op<T>::myNeg(a*v.m_g[i]);
	}
	else
	{
		USER_ASSERT(m_gsize==v.size(), ".diff(i,n) called with two different values of"
		" n in same expression; "<< m_gsize << " and " << v.size())
		for(int i=0;i<m_gsize;++i) Op<T>::myCsub(m_g[i],a*v.m_g[i]);
	}
}

// SPARSE BACKWARD NODE

template <class T>
INLINE2 BackwardNodeSparse<T>::~BackwardNodeSparse()
{
	ASSERT(m_gsize?m_g!=0:m_g==0)
	if (m_g)
	{
		for(int i=0;i<m_gsize;++i) delScalar(T,m_g[i]);
		delArray(T*,m_g);
	}
	m_g=0;
}
template <class T>
INLINE2 void BackwardNodeSparse<T>::touchg(int n)
{
	if (!m_g && n>0)
	{
		m_g=newArray(T*,n);
		m_gsize=n;
		for(int i=0;i<n;m_g[i++]=0);
	}
}
template <class T>
INLINE1 T& BackwardNodeSparse<T>::value(){ return m_v; }
template <class T>
INLINE1 T& BackwardNodeSparse<T>::deriv(const int idx){ return *(m_g[idx]?m_g[idx]:m_g[idx]=newCopy(T,0)); }
template <class T>
INLINE1 const T& BackwardNodeSparse<T>::deriv(const int idx) const { static const T z(0); return m_g[idx]?*m_g[idx]:z; }
template <class T>
INLINE1 int BackwardNodeSparse<T>::size() const { return m_gsize; }
template <class T>
INLINE1 int BackwardNodeSparse<T>::rc() const { return m_rc; }
template <class T>
INLINE1 int BackwardNodeSparse<T>::inc_rc() { return ++m_rc; }
template <class T>
INLINE1 int BackwardNodeSparse<T>::dec_rc() { return --m_rc; }
template <class T>
INLINE2 void BackwardNodeSparse<T>::pv(const BackwardNodeSparse& v)
{
	if (0==v.m_g) return;
	if (0==m_g)
	{
		m_gsize=v.size();
		m_g=newArray(T*,m_gsize);
		for(int i=0;i<m_gsize;m_g[i++]=0);
	}
	else
	{
		USER_ASSERT(m_gsize==v.size(), ".diff(i,n) called with two different values of"
		" n in same expression; "<< m_gsize << " and " << v.size())
	}
	for(int i=0;i<m_gsize;++i)
	{
		if (v.m_g[i])
		{
			if (m_g[i]) Op<T>::myCadd(*m_g[i],*v.m_g[i]);
			else m_g[i]=newCopy(T,*v.m_g[i]);
		}
	}
}
template <class T>
INLINE2 void BackwardNodeSparse<T>::mv(const BackwardNodeSparse& v)
{
	if (0==v.m_g) return;
	if (0==m_g)
	{
		m_gsize=v.size();
		m_g=newArray(T*,m_gsize);
		for(int i=0;i<m_gsize;m_g[i++]=0);
	}
	else
	{
		USER_ASSERT(m_gsize==v.size(), ".diff(i,n) called with two different values of"
		" n in same expression; "<< m_gsize << " and " << v.size())
	}
	for(int i=0;i<m_gsize;++i)
	{
		if (v.m_g[i])
		{
			if (m_g[i]) Op<T>::myCsub(*m_g[i], *v.m_g[i]);
			else m_g[i]=newCopy(T,Op<T>::myNeg(*v.m_g[i]));
		}
	}
}
template <class T>
INLINE2 void BackwardNodeSparse<T>::pav(const T& a, const BackwardNodeSparse& v)
{
	if (0==v.m_g) return;
	if (0==m_g)
	{
		m_gsize=v.size();
		m_g=newArray(T*,m_gsize);
		for(int i=0;i<m_gsize;m_g[i++]=0);
	}
	else
	{
		USER_ASSERT(m_gsize==v.size(), ".diff(i,n) called with two different values of"
		" n in same expression; "<< m_gsize << " and " << v.size())
	}
	for(int i=0;i<m_gsize;++i)
	{
		if (v.m_g[i])
		{
			if (m_g[i]) Op<T>::myCadd(*m_g[i],a**v.m_g[i]);
			else m_g[i]=newCopy(T,a**v.m_g[i]);
		}
	}
}
template <class T>
INLINE2 void BackwardNodeSparse<T>::mav(const T& a, const BackwardNodeSparse& v)
{
	if (0==v.m_g) return;
	if (0==m_g)
	{
		m_gsize=v.size();
		m_g=newArray(T*,m_gsize);
		for(int i=0;i<m_gsize;m_g[i++]=0);
	}
	else
	{
		USER_ASSERT(m_gsize==v.size(), ".diff(i,n) called with two different values of"
		" n in same expression; "<< m_gsize << " and " << v.size())
	}
	for(int i=0;i<m_gsize;++i)
	{
		if (v.m_g[i])
		{
			if (m_g[i]) Op<T>::myCsub(*m_g[i],a**v.m_g[i]);
			else m_g[i]=newCopy(T,Op<T>::myNeg(a**v.m_g[i]));
		}
	}
}

// BackwardOp

template <class T>
INLINE1 BackwardOp<T>::BackwardOp():m_val(0)
{
}

template <class T>
INLINE1 BackwardOp<T>::BackwardOp(const int x):m_val(newCopy(BackwardNode,x))
{
}

template <class T>
INLINE1 BackwardOp<T>::BackwardOp(const T &x):m_val(newCopy(BackwardNode,x))
{
}

template <class T>
INLINE1 BackwardOp<T>::BackwardOp(const BackwardOp<T> &x):m_val(x.m_val)
{
	if (m_val) inc_rc();
}

template <class T>
INLINE1 BackwardOp<T>::~BackwardOp()
{
	dec_rc();
}

template <class T>
INLINE1 void BackwardOp<T>::inc_rc() const
{
/* Increments the resource counter of the underlying BackwardNode. */
	ASSERT(m_val);
	m_val->inc_rc();
}

template <class T>
INLINE2 void BackwardOp<T>::dec_rc()
{
/* Decrements the resource counter of the underlying BackwardNode. */
	if (m_val && m_val->dec_rc()==0)
	{
	/* When no operator refers to a BackwardNode it can be deallocated. */
		delScalar(BackwardNode,m_val);
		m_val=0;
	}
}

/* Dummy functions */
template <class T>
INLINE1 BackwardOp<T>* BackwardOp<T>::copy() const
{
	return 0;
}

template <class T>
INLINE1 void BackwardOp<T>::propagate()
{
}

template <class T>
INLINE1 void BackwardOp<T>::propagateop()
{
}

template <class T>
INLINE1 bool operator == (const BackwardOp<T> &a, const BackwardOp<T> &b)
{
	ASSERT(a.m_val && b.m_val);
	return Op<T>::myEq(a.m_val->value(),b.m_val->value());
};

template <class T>
INLINE1 bool operator == (const BackwardOp<T> &a, const typename Op<T>::Base &b)
{
	ASSERT(a.m_val);
	return Op<T>::myEq(a.m_val->value(),b);
}

template <class T>
INLINE1 bool operator == (const typename Op<T>::Base &a, const BackwardOp<T> &b)
{
	ASSERT(b.m_val);
	return Op<T>::myEq(a,b.m_val->value());
};

template <class T>
INLINE1 bool operator != (const BackwardOp<T> &a, const BackwardOp<T> &b)
{
	ASSERT(a.m_val && b.m_val);
	return Op<T>::myNe(a.m_val->value(),b.m_val->value());
}

template <class T>
INLINE1 bool operator != (const BackwardOp<T> &a, const typename Op<T>::Base &b)
{
	ASSERT(a.m_val);
	return Op<T>::myNe(a.m_val->value(),b);
}

template <class T>
INLINE1 bool operator != (const typename Op<T>::Base &a, const BackwardOp<T> &b)
{
	ASSERT(b.m_val);
	return Op<T>::myNe(a,b.m_val->value());
}

template <class T>
INLINE1 bool operator > (const BackwardOp<T> &a, const BackwardOp<T> &b)
{
	ASSERT(a.m_val && b.m_val);
	return Op<T>::myGt(a.m_val->value(),b.m_val->value());
}

template <class T>
INLINE1 bool operator > (const BackwardOp<T> &a, const typename Op<T>::Base &b)
{
	ASSERT(a.m_val);
	return Op<T>::myGt(a.m_val->value(),b);
}

template <class T>
INLINE1 bool operator > (const typename Op<T>::Base &a, const BackwardOp<T> &b)
{
	ASSERT(b.m_val);
	return Op<T>::myGt(a,b.m_val->value());
}

template <class T>
INLINE1 bool operator >= (const BackwardOp<T> &a, const BackwardOp<T> &b)
{
	ASSERT(a.m_val && b.m_val);
	return Op<T>::myGe(a.m_val->value(),b.m_val->value());
}

template <class T>
INLINE1 bool operator >= (const BackwardOp<T> &a, const typename Op<T>::Base &b)
{
	ASSERT(a.m_val);
	return Op<T>::myGe(a.m_val->value(),b);
}

template <class T>
INLINE1 bool operator >= (const typename Op<T>::Base &a, const BackwardOp<T> &b)
{
	ASSERT(b.m_val);
	return Op<T>::myGe(a,b.m_val->value());
}

template <class T>
INLINE1 bool operator < (const BackwardOp<T> &a, const BackwardOp<T> &b)
{
	ASSERT(a.m_val && b.m_val);
	return Op<T>::myLt(a.m_val->value(),b.m_val->value());
}

template <class T>
INLINE1 bool operator < (const BackwardOp<T> &a, const typename Op<T>::Base &b)
{
	ASSERT(a.m_val);
	return Op<T>::myLt(a.m_val->value(),b);
}

template <class T>
INLINE1 bool operator < (const typename Op<T>::Base &a, const BackwardOp<T> &b)
{
	ASSERT(b.m_val);
	return Op<T>::myLt(a,b.m_val->value());
}

template <class T>
INLINE1 bool operator <= (const BackwardOp<T> &a, const BackwardOp<T> &b)
{
	ASSERT(a.m_val && b.m_val);
	return Op<T>::myLe(a.m_val->value(),b.m_val->value());
}

template <class T>
INLINE1 bool operator <= (const BackwardOp<T> &a, const typename Op<T>::Base &b)
{
	ASSERT(a.m_val);
	return Op<T>::myLe(a.m_val->value(),b);
}

template <class T>
INLINE1 bool operator <= (const typename Op<T>::Base &a, const BackwardOp<T> &b)
{
	ASSERT(b.m_val);
	return Op<T>::myLe(a,b.m_val->value());
}

template <class T>
INLINE1 BINBackwardOp<T>::BINBackwardOp():m_o1(0),m_o2(0)
{
}

template <class T>
INLINE1 BINBackwardOp<T>::BINBackwardOp(const T& x):BackwardOp<T>(x),m_o1(0),m_o2(0)
{
}

template <class T>
INLINE1 BINBackwardOp<T>::BINBackwardOp(const T& x,BackwardOp<T>* a,BackwardOp<T>* b):BackwardOp<T>(x),m_o1(a),m_o2(b)
{
}

template <class T>
INLINE1 BINBackwardOp<T>::BINBackwardOp(const BINBackwardOp<T>& x):BackwardOp<T>(x),m_o1(x.m_o1),m_o2(x.m_o2)
{
}

template <class T>
INLINE2 BINBackwardOp<T>::~BINBackwardOp()
{
	if (this->m_val && this->m_val->rc()==1)
	{
		if (m_o1) {m_o1->propagate(); delScalar(BackwardOp<T>,m_o1); m_o1=0;}
		if (m_o2) {m_o2->propagate(); delScalar(BackwardOp<T>,m_o2); m_o2=0;}
	}
}

template <class T>
INLINE2 void BINBackwardOp<T>::propagate()
{
// This function propagates the partial derivatives. The propagation
// does not take place until all incomming partial derivatives have
// been received.
//
	if (this->m_val && this->m_val->rc()==1)
	{
		this->propagateop();
		if (this->m_o1) {this->m_o1->propagate(); delScalar(BackwardOp<T>,this->m_o1); this->m_o1=0;}
		if (this->m_o2) {this->m_o2->propagate(); delScalar(BackwardOp<T>,this->m_o2); this->m_o2=0;}
	}
	this->dec_rc();
}

template <class T>
INLINE1 UNBackwardOp<T>::UNBackwardOp():m_o1(0)
{
}

template <class T>
INLINE1 UNBackwardOp<T>::UNBackwardOp(const T& x):BackwardOp<T>(x),m_o1(0)
{
}

template <class T>
INLINE1 UNBackwardOp<T>::UNBackwardOp(const int x):BackwardOp<T>(x),m_o1(0)
{
}

template <class T>
INLINE1 UNBackwardOp<T>::UNBackwardOp(const T& x,BackwardOp<T>* a):BackwardOp<T>(x),m_o1(a)
{
}

template <class T>
INLINE1 UNBackwardOp<T>::UNBackwardOp(const UNBackwardOp<T>& x):BackwardOp<T>(x),m_o1(x.m_o1)
{
}

template <class T>
INLINE2 UNBackwardOp<T>::~UNBackwardOp()
{
	if (this->m_val && this->m_val->rc()==1)
	{
		if (m_o1)
		{
			m_o1->propagate();
			delScalar(BackwardOp<T>,m_o1);
			m_o1=0;
		}
	}
}

template <class T>
INLINE2 void UNBackwardOp<T>::propagate()
{
// This function propagates the partial derivatives. The propagation
// does not take place until all incomming partial derivatives have
// been received.
//
	if (this->m_val)
	{
		ASSERT(this->m_val->rc()==1);
		this->propagateop();
		if (m_o1)
		{
			m_o1->propagate();
			delScalar(BackwardOp<T>,m_o1);
			m_o1=0;
		}
	}
	this->dec_rc();
}

template <class T> template <class U, class V>
BTypeName<T>::Assign<U,V>::Assign(BTypeName<U>& dest, const U& x)
{
	dest.dec_rc();
	dest.m_val=newCopy(typename BackwardOp<U>::BackwardNode,x);
	dest.m_o1=0;
}
template <class T> template <class U, class V>
BTypeName<T>::Assign<U,V>::Assign(BTypeName<U>& dest, const BTypeName<U>& x)
{
	dest.dec_rc();
	dest.m_val=x.m_val;
	if (dest.m_val) dest.inc_rc();
	dest.m_o1=x.m_o1;
}
template <class T> template <class U, class V>
BTypeName<T>::Assign<U,V>::Assign(BTypeName<U>& dest, const BackwardOp<U>& x)
{
	dest.dec_rc();
	dest.m_val=x.m_val;
	dest.inc_rc();
	dest.m_o1=x.copy();
}
template <class T> template <class U, class V>
BTypeName<T>::Assign<U,V>::Assign(BTypeName<U>& dest, const typename Op<U>::Base& x)
{
	dest.dec_rc();
	dest.m_val=newCopy(typename BackwardOp<T>::BackwardNode,x);
	dest.m_o1=0;
}

template <class T> template <class U>
BTypeName<T>::Assign<U,U>::Assign(BTypeName<U>& dest, const U& x)
{
	dest.dec_rc();
	dest.m_val=newCopy(typename BackwardOp<U>::BackwardNode,x);
	dest.m_o1=0;
}
template <class T> template <class U>
BTypeName<T>::Assign<U,U>::Assign(BTypeName<U>& dest, const BTypeName<U>& x)
{
	dest.dec_rc();
	dest.m_val=x.m_val;
	if (dest.m_val) dest.inc_rc();
	dest.m_o1=x.m_o1;
}
template <class T> template <class U>
BTypeName<T>::Assign<U,U>::Assign(BTypeName<U>& dest, const BackwardOp<U>& x)
{
	dest.dec_rc();
	dest.m_val=x.m_val;
	dest.inc_rc();
	dest.m_o1=x.copy();
}

template <class T>
INLINE1 BTypeName<T>::BTypeName()
{
}

template <class T>
INLINE1 BTypeName<T>::BTypeName(const int x):UNBackwardOp<T>(x)
{
}

template <class T>
INLINE1 BTypeName<T>::BTypeName(const BTypeName<T>& x):UNBackwardOp<T>(x)
{
}

template <class T>
INLINE2 BTypeName<T>::BTypeName(const BackwardOp<T>& x)
{
	this->m_val=x.m_val;
	this->m_o1=x.copy();
	this->inc_rc();
}

template <class T> template <class U> 
INLINE2 BTypeName<T>::BTypeName(const U& x)
{
	Assign< T, typename Op<T>::Base >(*this,x);
}



template <class T>
INLINE2 BTypeName<T>::~BTypeName()
{
	if (this->m_val && this->m_val->rc()==2 && this->m_o1)
	{
	// If object is the last BTypeName encapsulating some operator then
	// underlying operator should also be deleted.
	//
		this->dec_rc(); // When called again in ~op, BackwardNode is already deleted.
	}
}

template <class T>
INLINE2 void BTypeName<T>::dec_rc()
{
// If this is the last BTypeName encapsulating some operator then operator
// can be deleted (and this will delete the operators BackwardNode).
//
	if (this->m_o1)
	{
		ASSERT(this->m_val);
		if (this->m_val->dec_rc()==1)
		{
			delScalar(BackwardOp<T>,this->m_o1);
			this->m_o1=0;
			this->m_val=0; // Since deleted from o1.
		}
	} else UNBackwardOp<T>::dec_rc();
}

template <class T>
INLINE2 void BTypeName<T>::propagate()
{
// This function propagates the partial derivatives. The propagation
// does not take place until all incomming partial derivatives have
// been received. If object is the last object encapsulating some
// operator then propagation has to be passed on to this.
//
	ASSERT( this->m_val || this->m_o1==0 );
	if (this->m_o1 && this->m_val->rc()==2)
	{
		UNBackwardOp<T>::dec_rc();
		this->m_o1->propagate();
		delScalar(BackwardOp<T>,this->m_o1);
	} else
		UNBackwardOp<T>::dec_rc();
	this->m_o1=0;
	this->m_val=0;
}

template <class T>
INLINE2 void BTypeName<T>::diff(int idx, int n)
{
// Initiates the propagation of derivatives, ie. the backward sweep.
// idx is the index of the variable this objects represent of all
// the n "dependent" variables.
//
	ASSERT(this->m_val);
	this->m_val->touchg(n);
	this->m_val->deriv(idx)=Op<T>::myOne();
	BTypeName<T> newthis(this->m_val->value());
	this->propagate();
	(*this)=newthis;
	this->m_o1 = 0;
}

template <class T>
INLINE2 BackwardOp<T>* BTypeName<T>::copy() const
{
// Returns a copy of the object which still links to the same copy
// of the internal BackwardNode.
//
	return newCopy(BTypeName<T>,*this);
}

template <class T>
INLINE2 T& BTypeName<T>::x()
{
	if (this->m_val)
		return this->m_val->value();
	else
	{
		static T zero(Op<T>::myZero());
		return zero;
	}
}

template <class T>
INLINE2 T& BTypeName<T>::d(int i)
{
	if (this->m_val && this->m_val->size())
	{
		ASSERT(i<this->m_val->size());
		return this->m_val->deriv(i);
	}
	else
	{
		static T zero(Op<T>::myZero());
		return zero;
	}
}

/* --------------------------------------------------------------------- */
/*                                                                       */
/*              BACKWARD DIFFERENTIATION ADDITION FORMULAS               */
/*                                                                       */
/* --------------------------------------------------------------------- */

/* --------------------------- op + op --------------------------------- */

template <class T>
INLINE1 ADDBackwardOp<T>::ADDBackwardOp(const ADDBackwardOp<T>& o):
	BINBackwardOp<T>(o)
{}

template <class T>
INLINE1 ADDBackwardOp<T>::ADDBackwardOp(const BackwardOp<T>& a, const BackwardOp<T>& b):
	BINBackwardOp<T>(a.m_val->value()+b.m_val->value(),a.copy(),b.copy())
{}

template <class T>
INLINE1 ADDBackwardOp<T>::~ADDBackwardOp()
{}

template <class T>
INLINE2 BackwardOp<T>* ADDBackwardOp<T>::copy() const
{
	return newCopy(ADDBackwardOp<T>,*this);
}

template <class T>
INLINE2 void ADDBackwardOp<T>::propagateop()
{
	this->m_o1->m_val->pv(*this->m_val);
	this->m_o2->m_val->pv(*this->m_val);
}

/* --------------------------- const + op ------------------------------ */

template <class T, class U>
INLINE1 ADDBackwardOp3<T,U>::ADDBackwardOp3(const ADDBackwardOp3<T,U>& o):
	UNBackwardOp<T>(o),m_ax(o.m_ax)
{}

template <class T, class U>
INLINE1 ADDBackwardOp3<T,U>::ADDBackwardOp3(const U& a, const BackwardOp<T>& b):
	UNBackwardOp<T>(a+b.m_val->value())
{
	this->m_o1=b.copy();
	m_ax=a;
}

template <class T, class U>
INLINE1 ADDBackwardOp3<T,U>::~ADDBackwardOp3()
{}

template <class T, class U>
INLINE2 BackwardOp<T>* ADDBackwardOp3<T,U>::copy() const
{
	typedef ADDBackwardOp3<T,U> MyType;
	return newCopy(MyType,*this);
}

template <class T, class U>
INLINE2 void ADDBackwardOp3<T,U>::propagateop()
{
	this->m_o1->m_val->pv(*this->m_val);
}

/* --------------------------- op + const ------------------------------ */

template <class T, class U>
INLINE1 ADDBackwardOp4<T,U>::ADDBackwardOp4(const ADDBackwardOp4<T,U>& o):
	UNBackwardOp<T>(o),m_bx(o.m_bx)
{}

template <class T, class U>
INLINE1 ADDBackwardOp4<T,U>::ADDBackwardOp4(const BackwardOp<T>& a, const U& b):
	UNBackwardOp<T>(a.m_val->value()+b)
{
	this->m_o1=a.copy();
	m_bx=b;
}

template <class T, class U>
INLINE1 ADDBackwardOp4<T,U>::~ADDBackwardOp4()
{}

template <class T, class U>
INLINE2 BackwardOp<T>* ADDBackwardOp4<T,U>::copy() const
{
	typedef ADDBackwardOp4<T,U> MyType;
	return newCopy(MyType,*this);
}

template <class T, class U>
INLINE2 void ADDBackwardOp4<T,U>::propagateop()
{
	this->m_o1->m_val->pv(*this->m_val);
}

/* --------------------------------------------------------------------- */
/*                                                                       */
/*            BACKWARD DIFFERENTIATION SUBTRACTION FORMULAS              */
/*                                                                       */
/* --------------------------------------------------------------------- */

/* --------------------------- op - op --------------------------------- */

template <class T>
INLINE1 SUBBackwardOp<T>::SUBBackwardOp(const SUBBackwardOp<T>& o):
	BINBackwardOp<T>(o)
{}

template <class T>
INLINE1 SUBBackwardOp<T>::SUBBackwardOp(const BackwardOp<T>& a, const BackwardOp<T>& b):
	BINBackwardOp<T>(a.m_val->value()-b.m_val->value(),a.copy(),b.copy())
{}

template <class T>
INLINE1 SUBBackwardOp<T>::~SUBBackwardOp()
{}

template <class T>
INLINE2 BackwardOp<T>* SUBBackwardOp<T>::copy() const
{
	return newCopy(SUBBackwardOp<T>,*this);
}

template <class T>
INLINE2 void SUBBackwardOp<T>::propagateop()
{
	this->m_o1->m_val->pv(*this->m_val);
	this->m_o2->m_val->mv(*this->m_val);
}

/* --------------------------- const - op ------------------------------ */

template <class T, class U>
INLINE1 SUBBackwardOp3<T,U>::SUBBackwardOp3(const SUBBackwardOp3<T,U>& o):
	UNBackwardOp<T>(o),m_ax(o.m_ax)
{}

template <class T, class U>
INLINE1 SUBBackwardOp3<T,U>::SUBBackwardOp3(const U& a, const BackwardOp<T>& b):
	UNBackwardOp<T>(a-b.m_val->value())
{
	this->m_o1=b.copy();
	m_ax=a;
}

template <class T, class U>
INLINE1 SUBBackwardOp3<T,U>::~SUBBackwardOp3()
{}

template <class T, class U>
INLINE2 BackwardOp<T>* SUBBackwardOp3<T,U>::copy() const
{
	typedef SUBBackwardOp3<T,U> MyType;
	return newCopy(MyType,*this);
}

template <class T, class U>
INLINE2 void SUBBackwardOp3<T,U>::propagateop()
{
	this->m_o1->m_val->mv(*this->m_val);
}

/* --------------------------- op - const ------------------------------ */

template <class T, class U>
INLINE1 SUBBackwardOp4<T,U>::SUBBackwardOp4(const SUBBackwardOp4<T,U>& o):
	UNBackwardOp<T>(o),m_bx(o.m_bx)
{}

template <class T, class U>
INLINE1 SUBBackwardOp4<T,U>::SUBBackwardOp4(const BackwardOp<T>& a, const U& b):
	UNBackwardOp<T>(a.m_val->value()-b)
{
	this->m_o1=a.copy();
	m_bx=b;
}

template <class T, class U>
INLINE1 SUBBackwardOp4<T,U>::~SUBBackwardOp4()
{}

template <class T, class U>
INLINE2 BackwardOp<T>* SUBBackwardOp4<T,U>::copy() const
{
	typedef SUBBackwardOp4<T,U> MyType;
	return newCopy(MyType,*this);
}

template <class T, class U>
INLINE2 void SUBBackwardOp4<T,U>::propagateop()
{
	this->m_o1->m_val->pv(*this->m_val);
}

/* --------------------------------------------------------------------- */
/*                                                                       */
/*          BACKWARD DIFFERENTIATION MULTIPLICATION FORMULAS             */
/*                                                                       */
/* --------------------------------------------------------------------- */

/* --------------------------- op * op --------------------------------- */

template <class T>
INLINE1 MULBackwardOp<T>::MULBackwardOp(const MULBackwardOp<T>& o):
	BINBackwardOp<T>(o)
{}

template <class T>
INLINE1 MULBackwardOp<T>::MULBackwardOp(const BackwardOp<T>& a, const BackwardOp<T>& b):
	BINBackwardOp<T>(a.m_val->value()*b.m_val->value(),a.copy(),b.copy())
{}

template <class T>
INLINE1 MULBackwardOp<T>::~MULBackwardOp()
{}

template <class T>
INLINE2 BackwardOp<T>* MULBackwardOp<T>::copy() const
{
	return newCopy(MULBackwardOp<T>,*this);
}

template <class T>
INLINE2 void MULBackwardOp<T>::propagateop()
{
	this->m_o1->m_val->pav(this->m_o2->m_val->value(),*this->m_val);
	this->m_o2->m_val->pav(this->m_o1->m_val->value(),*this->m_val);
}

/* --------------------------- const * op ------------------------------ */

template <class T, class U>
INLINE1 MULBackwardOp3<T,U>::MULBackwardOp3(const MULBackwardOp3<T,U>& o):
	UNBackwardOp<T>(o),m_ax(o.m_ax)
{}

template <class T, class U>
INLINE1 MULBackwardOp3<T,U>::MULBackwardOp3(const U& a, const BackwardOp<T>& b):
	UNBackwardOp<T>(a*b.m_val->value())
{
	this->m_o1=b.copy();
	this->m_ax=a;
}

template <class T, class U>
INLINE1 MULBackwardOp3<T,U>::~MULBackwardOp3()
{}

template <class T, class U>
INLINE2 BackwardOp<T>* MULBackwardOp3<T,U>::copy() const
{
	typedef MULBackwardOp3<T,U> MyType;
	return newCopy(MyType,*this);
}

template <class T, class U>
INLINE2 void MULBackwardOp3<T,U>::propagateop()
{
	this->m_o1->m_val->pav(T(m_ax),*this->m_val);
}

/* --------------------------- op * const ------------------------------ */

template <class T, class U>
INLINE1 MULBackwardOp4<T,U>::MULBackwardOp4(const MULBackwardOp4<T,U>& o):
	UNBackwardOp<T>(o),m_bx(o.m_bx)
{}

template <class T, class U>
INLINE1 MULBackwardOp4<T,U>::MULBackwardOp4(const BackwardOp<T>& a, const U& b):
	UNBackwardOp<T>(a.m_val->value()*b)
{
	this->m_o1=a.copy();
	this->m_bx=b;
}

template <class T, class U>
INLINE1 MULBackwardOp4<T,U>::~MULBackwardOp4()
{}

template <class T, class U>
INLINE2 BackwardOp<T>* MULBackwardOp4<T,U>::copy() const
{
	typedef MULBackwardOp4<T,U> MyType;
	return newCopy(MyType,*this);
}

template <class T, class U>
INLINE2 void MULBackwardOp4<T,U>::propagateop()
{
	this->m_o1->m_val->pav(T(m_bx),*this->m_val);
}

/* --------------------------------------------------------------------- */
/*                                                                       */
/*              BACKWARD DIFFERENTIATION DIVISION FORMULAS               */
/*                                                                       */
/* --------------------------------------------------------------------- */

/* --------------------------- op / op --------------------------------- */

template <class T>
INLINE1 DIVBackwardOp<T>::DIVBackwardOp(const DIVBackwardOp<T>& o):
	BINBackwardOp<T>(o)
{}

template <class T>
INLINE1 DIVBackwardOp<T>::DIVBackwardOp(const BackwardOp<T>& a, const BackwardOp<T>& b):
	BINBackwardOp<T>(a.m_val->value()/b.m_val->value(),a.copy(),b.copy())
{}

template <class T>
INLINE1 DIVBackwardOp<T>::~DIVBackwardOp()
{}

template <class T>
INLINE2 BackwardOp<T>* DIVBackwardOp<T>::copy() const
{
	return newCopy(DIVBackwardOp<T>,*this);
}

template <class T>
INLINE2 void DIVBackwardOp<T>::propagateop()
{
	T tmp=Op<T>::myInv(this->m_o2->m_val->value());
	this->m_o1->m_val->pav(tmp,*this->m_val);
	this->m_o2->m_val->mav(tmp*this->m_val->value(),*this->m_val);
}

/* --------------------------- typename Op<T>::Base / op ------------------------------ */

template <class T, class U>
INLINE1 DIVBackwardOp3<T,U>::DIVBackwardOp3(const DIVBackwardOp3<T,U>& o):
	UNBackwardOp<T>(o),m_ax(o.m_ax)
{}

template <class T, class U>
INLINE1 DIVBackwardOp3<T,U>::DIVBackwardOp3(const U& a, const BackwardOp<T>& b):
	UNBackwardOp<T>(a/b.m_val->value())
{
	this->m_o1=b.copy();
	m_ax=a;
}

template <class T, class U>
INLINE1 DIVBackwardOp3<T,U>::~DIVBackwardOp3()
{}

template <class T, class U>
INLINE2 BackwardOp<T>* DIVBackwardOp3<T,U>::copy() const
{
	typedef DIVBackwardOp3<T,U> MyType;
	return newCopy(MyType,*this);
}

template <class T, class U>
INLINE2 void DIVBackwardOp3<T,U>::propagateop()
{
	T tmp(this->m_val->value() / (this->m_o1->m_val)->value());
	this->m_o1->m_val->mav(tmp,*this->m_val);
}

/* --------------------------- op / const ------------------------------ */

template <class T, class U>
INLINE1 DIVBackwardOp4<T,U>::DIVBackwardOp4(const DIVBackwardOp4& o):
	UNBackwardOp<T>(o),m_bx(o.m_bx)
{}

template <class T, class U>
INLINE1 DIVBackwardOp4<T,U>::DIVBackwardOp4(const BackwardOp<T>& a, const U& b):
	UNBackwardOp<T>(a.m_val->value()/b)
{
	this->m_o1=a.copy();
	m_bx=b;
}

template <class T, class U>
INLINE1 DIVBackwardOp4<T,U>::~DIVBackwardOp4()
{}

template <class T, class U>
INLINE2 BackwardOp<T>* DIVBackwardOp4<T,U>::copy() const
{
	typedef DIVBackwardOp4<T,U> MyType;
	return newCopy(MyType,*this);
}

template <class T, class U>
INLINE2 void DIVBackwardOp4<T,U>::propagateop()
{
	this->m_o1->m_val->pav(T(Op<U>::myInv(m_bx)),*this->m_val);
}

/* --------------------------------------------------------------------- */
/*                                                                       */
/*                BACKWARD DIFFERENTIATION POWER FORMULAS                */
/*                                                                       */
/* --------------------------------------------------------------------- */

/* -------------------------- pow(op,BackwardOp) ------------------------------- */

template <class T>
INLINE1 POWBackwardOp<T>::POWBackwardOp(const POWBackwardOp& o):
	BINBackwardOp<T>(o)
{}

template <class T>
INLINE1 POWBackwardOp<T>::POWBackwardOp(const BackwardOp<T>& a, const BackwardOp<T>& b):
	BINBackwardOp<T>(Op<T>::myPow(a.m_val->value(),b.m_val->value()),a.copy(),b.copy())
{}

template <class T>
INLINE1 POWBackwardOp<T>::~POWBackwardOp()
{}

template <class T>
INLINE2 BackwardOp<T>* POWBackwardOp<T>::copy() const
{
	return newCopy(POWBackwardOp<T>,*this);
}

template <class T>
INLINE2 void POWBackwardOp<T>::propagateop()
{
	T tmp1(this->m_o2->m_val->value() *  Op<T>::myPow(this->m_o1->m_val->value(),this->m_o2->m_val->value()-Op<T>::myOne()));
	T tmp2(this->m_val->value() * Op<T>::myLog(this->m_o1->m_val->value()));
	this->m_o1->m_val->pav(tmp1,*this->m_val);
	this->m_o2->m_val->pav(tmp2,*this->m_val);
}

/* ----------------------- pow(typename Op<T>::Base,BackwardOp) ---------------------------- */

template <class T>
INLINE1 POWBackwardOp3<T>::POWBackwardOp3(const POWBackwardOp3& o):
	UNBackwardOp<T>(o)
{}

template <class T>
INLINE1 POWBackwardOp3<T>::POWBackwardOp3(const typename Op<T>::Base& a, const BackwardOp<T>& b):
	UNBackwardOp<T>(Op<T>::myPow(a,b.m_val->value()))
{
	this->m_o1=b.copy();
	m_ax=a;
}

template <class T>
INLINE1 POWBackwardOp3<T>::~POWBackwardOp3()
{}

template <class T>
INLINE2 BackwardOp<T>* POWBackwardOp3<T>::copy() const
{
	return newCopy(POWBackwardOp3<T>,*this);
}

template <class T>
INLINE2 void POWBackwardOp3<T>::propagateop()
{
	T tmp(this->m_val->value() * Op<T>::myLog(m_ax));
	this->m_o1->m_val->pav(tmp,*this->m_val);
}

/* ------------------------ pow(op,typename Op<T>::Base) --------------------------- */

template <class T>
INLINE1 POWBackwardOp4<T>::POWBackwardOp4(const POWBackwardOp4& o):
	UNBackwardOp<T>(o),m_bx(o.m_bx)
{}

template <class T>
INLINE1 POWBackwardOp4<T>::POWBackwardOp4(const BackwardOp<T>& a, const typename Op<T>::Base& b):
	UNBackwardOp<T>(Op<T>::myPow(a.m_val->value(),b))
{
	this->m_o1=a.copy();
	m_bx=b;
}

template <class T>
INLINE1 POWBackwardOp4<T>::~POWBackwardOp4()
{}

template <class T>
INLINE2 BackwardOp<T>* POWBackwardOp4<T>::copy() const
{
	return newCopy(POWBackwardOp4<T>,*this);
}

template <class T>
INLINE2 void POWBackwardOp4<T>::propagateop()
{
	T tmp(m_bx * Op<T>::myPow(this->m_o1->m_val->value(),m_bx-Op<T>::myOne()));
	this->m_o1->m_val->pav(tmp,*this->m_val);
}

/* --------------------------------------------------------------------- */
/*                                                                       */
/*            BACKWARD DIFFERENTIATION UNARY + FORMULA                   */
/*                                                                       */
/* --------------------------------------------------------------------- */

template <class T>
INLINE1 UADDBackwardOp<T>::UADDBackwardOp(const UADDBackwardOp& o):
	UNBackwardOp<T>(o)
{}

template <class T>
INLINE1 UADDBackwardOp<T>::UADDBackwardOp(const BackwardOp<T>& a):
	UNBackwardOp<T>(a.m_val->value(),a.copy())
{}

template <class T>
INLINE1 UADDBackwardOp<T>::~UADDBackwardOp()
{}

template <class T>
INLINE2 BackwardOp<T>* UADDBackwardOp<T>::copy() const
{
	return newCopy(UADDBackwardOp<T>,*this);
}

template <class T>
INLINE2 void UADDBackwardOp<T>::propagateop()
{
	this->m_o1->m_val->pv(*this->m_val);
}

/* --------------------------------------------------------------------- */
/*                                                                       */
/*            BACKWARD DIFFERENTIATION UNARY - FORMULA                   */
/*                                                                       */
/* --------------------------------------------------------------------- */

template <class T>
INLINE1 USUBBackwardOp<T>::USUBBackwardOp(const USUBBackwardOp& o):
	UNBackwardOp<T>(o)
{}

template <class T>
INLINE1 USUBBackwardOp<T>::USUBBackwardOp(const BackwardOp<T>& a):
	UNBackwardOp<T>(Op<T>::myNeg(a.m_val->value()),a.copy())
{}

template <class T>
INLINE1 USUBBackwardOp<T>::~USUBBackwardOp()
{}

template <class T>
INLINE2 BackwardOp<T>* USUBBackwardOp<T>::copy() const
{
	return newCopy(USUBBackwardOp<T>,*this);
}

template <class T>
INLINE2 void USUBBackwardOp<T>::propagateop()
{
	this->m_o1->m_val->mv(*this->m_val);
}

/* --------------------------------------------------------------------- */
/*                                                                       */
/*          BACKWARD DIFFERENTIATION pow(op,int) FORMULA                 */
/*                                                                       */
/* --------------------------------------------------------------------- */

template <class T>
INLINE1 POWopN<T>::POWopN(const POWopN& o):
	UNBackwardOp<T>(o),m_b(o.m_b)
{}

template <class T>
INLINE1 POWopN<T>::POWopN(const BackwardOp<T>& a, int ib):
	UNBackwardOp<T>(Op<T>::myPow(a.m_val->value(),ib),a.copy()),m_b(ib)
{}

template <class T>
INLINE1 POWopN<T>::~POWopN()
{}

template <class T>
INLINE2 BackwardOp<T>* POWopN<T>::copy() const
{
	return newCopy(POWopN<T>,*this);
}

template <class T>
INLINE2 void POWopN<T>::propagateop()
{
	T tmp(m_b*Op<T>::myPow(this->m_o1->m_val->value(),m_b-1));
	this->m_o1->m_val->pav(tmp,*this->m_val);
}

/* --------------------------------------------------------------------- */
/*                                                                       */
/*              BACKWARD DIFFERENTIATION sqr FORMULA                     */
/*                                                                       */
/* --------------------------------------------------------------------- */

template <class T>
INLINE1 SQRBackwardOp<T>::SQRBackwardOp(const SQRBackwardOp& o):
	UNBackwardOp<T>(o)
{}

template <class T>
INLINE1 SQRBackwardOp<T>::SQRBackwardOp(const BackwardOp<T>& a):
	UNBackwardOp<T>(Op<T>::mySqr(a.m_val->value()),a.copy())
{}

template <class T>
INLINE1 SQRBackwardOp<T>::~SQRBackwardOp()
{}

template <class T>
INLINE2 BackwardOp<T>* SQRBackwardOp<T>::copy() const
{
	return newCopy(SQRBackwardOp<T>,*this);
}

template <class T>
INLINE2 void SQRBackwardOp<T>::propagateop()
{
	T tmp(Op<T>::myTwo()*(this->m_o1->m_val)->value());
	this->m_o1->m_val->pav(tmp,*this->m_val);
}

/* --------------------------------------------------------------------- */
/*                                                                       */
/*              BACKWARD DIFFERENTIATION exp FORMULA                     */
/*                                                                       */
/* --------------------------------------------------------------------- */

template <class T>
INLINE1 EXPBackwardOp<T>::EXPBackwardOp(const EXPBackwardOp& o):
	UNBackwardOp<T>(o)
{}

template <class T>
INLINE1 EXPBackwardOp<T>::EXPBackwardOp(const BackwardOp<T>& a):
	UNBackwardOp<T>(Op<T>::myExp(a.m_val->value()),a.copy())
{}

template <class T>
INLINE1 EXPBackwardOp<T>::~EXPBackwardOp()
{}

template <class T>
BackwardOp<T>* EXPBackwardOp<T>::copy() const
{
	return newCopy(EXPBackwardOp<T>,*this);
}

template <class T>
void EXPBackwardOp<T>::propagateop()
{
	this->m_o1->m_val->pav(this->m_val->value(),*this->m_val);
}

/* --------------------------------------------------------------------- */
/*                                                                       */
/*              BACKWARD DIFFERENTIATION log FORMULA                     */
/*                                                                       */
/* --------------------------------------------------------------------- */

template <class T>
INLINE1 LOGBackwardOp<T>::LOGBackwardOp(const LOGBackwardOp& o):
	UNBackwardOp<T>(o)
{}

template <class T>
INLINE1 LOGBackwardOp<T>::LOGBackwardOp(const BackwardOp<T>& a):
	UNBackwardOp<T>(Op<T>::myLog(a.m_val->value()),a.copy())
{}

template <class T>
INLINE1 LOGBackwardOp<T>::~LOGBackwardOp()
{}

template <class T>
BackwardOp<T>* LOGBackwardOp<T>::copy() const
{
	return newCopy(LOGBackwardOp<T>,*this);
}

template <class T>
void LOGBackwardOp<T>::propagateop()
{
	this->m_o1->m_val->pav(Op<T>::myInv(this->m_o1->m_val->value()),*this->m_val);
}

/* --------------------------------------------------------------------- */
/*                                                                       */
/*              BACKWARD DIFFERENTIATION sqrt FORMULA                    */
/*                                                                       */
/* --------------------------------------------------------------------- */

template <class T>
INLINE1 SQRTBackwardOp<T>::SQRTBackwardOp(const SQRTBackwardOp& o):
	UNBackwardOp<T>(o)
{}

template <class T>
INLINE1 SQRTBackwardOp<T>::SQRTBackwardOp(const BackwardOp<T>& a):
	UNBackwardOp<T>(Op<T>::mySqrt(a.m_val->value()),a.copy())
{}

template <class T>
INLINE1 SQRTBackwardOp<T>::~SQRTBackwardOp()
{}

template <class T>
INLINE2 BackwardOp<T>* SQRTBackwardOp<T>::copy() const
{
	return newCopy(SQRTBackwardOp<T>,*this);
}

template <class T>
INLINE2 void SQRTBackwardOp<T>::propagateop()
{
	T tmp(Op<T>::myInv(this->m_val->value()*Op<T>::myTwo()));
	this->m_o1->m_val->pav(tmp,*this->m_val);
}

/* --------------------------------------------------------------------- */
/*                                                                       */
/*              BACKWARD DIFFERENTIATION sin FORMULA                     */
/*                                                                       */
/* --------------------------------------------------------------------- */

template <class T>
INLINE1 SINBackwardOp<T>::SINBackwardOp(const SINBackwardOp& o):
	UNBackwardOp<T>(o)
{}

template <class T>
INLINE1 SINBackwardOp<T>::SINBackwardOp(const BackwardOp<T>& a):
	UNBackwardOp<T>(Op<T>::mySin(a.m_val->value()),a.copy())
{}

template <class T>
INLINE1 SINBackwardOp<T>::~SINBackwardOp()
{}

template <class T>
INLINE2 BackwardOp<T>* SINBackwardOp<T>::copy() const
{
	return newCopy(SINBackwardOp<T>,*this);
}

template <class T>
INLINE2 void SINBackwardOp<T>::propagateop()
{
	T tmp(Op<T>::myCos((this->m_o1->m_val)->value()));
	this->m_o1->m_val->pav(tmp,*this->m_val);
}

/* --------------------------------------------------------------------- */
/*                                                                       */
/*              BACKWARD DIFFERENTIATION cos FORMULA                     */
/*                                                                       */
/* --------------------------------------------------------------------- */

template <class T>
INLINE1 COSBackwardOp<T>::COSBackwardOp(const COSBackwardOp& o):
	UNBackwardOp<T>(o)
{}

template <class T>
INLINE1 COSBackwardOp<T>::COSBackwardOp(const BackwardOp<T>& a):
	UNBackwardOp<T>(Op<T>::myCos(a.m_val->value()),a.copy())
{}

template <class T>
INLINE1 COSBackwardOp<T>::~COSBackwardOp()
{}

template <class T>
INLINE2 BackwardOp<T>* COSBackwardOp<T>::copy() const
{
	return newCopy(COSBackwardOp<T>,*this);
}

template <class T>
INLINE2 void COSBackwardOp<T>::propagateop()
{
	T tmp(Op<T>::mySin((this->m_o1->m_val)->value()));
	this->m_o1->m_val->mav(tmp,*this->m_val);
}

/* --------------------------------------------------------------------- */
/*                                                                       */
/*              BACKWARD DIFFERENTIATION tan FORMULA                     */
/*                                                                       */
/* --------------------------------------------------------------------- */

template <class T>
INLINE1 TANBackwardOp<T>::TANBackwardOp(const TANBackwardOp& o):
	UNBackwardOp<T>(o)
{}

template <class T>
INLINE1 TANBackwardOp<T>::TANBackwardOp(const BackwardOp<T>& a):
	UNBackwardOp<T>(Op<T>::myTan(a.m_val->value()),a.copy())
{}

template <class T>
INLINE1 TANBackwardOp<T>::~TANBackwardOp()
{}

template <class T>
INLINE2 BackwardOp<T>* TANBackwardOp<T>::copy() const
{
	return newCopy(TANBackwardOp<T>,*this);
}

template <class T>
INLINE2 void TANBackwardOp<T>::propagateop()
{
	T tmp(Op<T>::mySqr(this->m_val->value())+Op<T>::myOne());
	this->m_o1->m_val->pav(tmp,*this->m_val);
}

/* --------------------------------------------------------------------- */
/*                                                                       */
/*              BACKWARD DIFFERENTIATION asin FORMULA                    */
/*                                                                       */
/* --------------------------------------------------------------------- */

template <class T>
INLINE1 ASINBackwardOp<T>::ASINBackwardOp(const ASINBackwardOp& o):
	UNBackwardOp<T>(o)
{}

template <class T>
INLINE1 ASINBackwardOp<T>::ASINBackwardOp(const BackwardOp<T>& a):
	UNBackwardOp<T>(Op<T>::myAsin(a.m_val->value()),a.copy())
{}

template <class T>
INLINE1 ASINBackwardOp<T>::~ASINBackwardOp()
{}

template <class T>
INLINE2 BackwardOp<T>* ASINBackwardOp<T>::copy() const
{
	return newCopy(ASINBackwardOp<T>,*this);
}

template <class T>
INLINE2 void ASINBackwardOp<T>::propagateop()
{
	T tmp(Op<T>::myInv(Op<T>::mySqrt(Op<T>::myOne()-Op<T>::mySqr(this->m_o1->m_val->value()))));
	this->m_o1->m_val->pav(tmp,*this->m_val);
}

/* --------------------------------------------------------------------- */
/*                                                                       */
/*              BACKWARD DIFFERENTIATION acos FORMULA                    */
/*                                                                       */
/* --------------------------------------------------------------------- */

template <class T>
INLINE1 ACOSBackwardOp<T>::ACOSBackwardOp(const ACOSBackwardOp& o):
	UNBackwardOp<T>(o)
{}

template <class T>
INLINE1 ACOSBackwardOp<T>::ACOSBackwardOp(const BackwardOp<T>& a):
	UNBackwardOp<T>(Op<T>::myAcos(a.m_val->value()),a.copy())
{}

template <class T>
INLINE1 ACOSBackwardOp<T>::~ACOSBackwardOp()
{}

template <class T>
INLINE2 BackwardOp<T>* ACOSBackwardOp<T>::copy() const
{
	return newCopy(ACOSBackwardOp<T>,*this);
}

template <class T>
INLINE2 void ACOSBackwardOp<T>::propagateop()
{
	T tmp(Op<T>::myInv(Op<T>::mySqrt(Op<T>::myOne()-Op<T>::mySqr(this->m_o1->m_val->value()))));
	this->m_o1->m_val->mav(tmp,*this->m_val);
}

/* --------------------------------------------------------------------- */
/*                                                                       */
/*              BACKWARD DIFFERENTIATION atan FORMULA                    */
/*                                                                       */
/* --------------------------------------------------------------------- */

template <class T>
ATANBackwardOp<T>::ATANBackwardOp(const ATANBackwardOp& o):
	UNBackwardOp<T>(o)
{}

template <class T>
INLINE1 ATANBackwardOp<T>::ATANBackwardOp(const BackwardOp<T>& a):
	UNBackwardOp<T>(Op<T>::myAtan(a.m_val->value()),a.copy())
{}

template <class T>
INLINE1 ATANBackwardOp<T>::~ATANBackwardOp()
{}

template <class T>
INLINE2 BackwardOp<T>* ATANBackwardOp<T>::copy() const
{
	return newCopy(ATANBackwardOp<T>,*this);
}

template <class T>
INLINE2 void ATANBackwardOp<T>::propagateop()
{
	T tmp(Op<T>::myInv(Op<T>::mySqr(this->m_o1->m_val->value())+Op<T>::myOne()));
	this->m_o1->m_val->pav(tmp,*this->m_val);
}

/* --------------------------------------------------------------------- */
/*                                                                       */
/*                        DEFINITION OF OPERATORS                        */
/*                                                                       */
/* --------------------------------------------------------------------- */

template <class T>
INLINE1 ADDBackwardOp<T> operator+ (const BackwardOp<T>& a, const BackwardOp<T>& b)
{
	ASSERT(a.m_val && b.m_val);
	return ADDBackwardOp<T>(a,b);
}

template <class T>
INLINE1 ADDBackwardOp<T> operator+ (const BackwardOp<T>& a, const BTypeName<T>& b)
{
	ASSERT(a.m_val && b.m_val);
	return ADDBackwardOp<T>(a,b);
}

template <class T>
INLINE1 ADDBackwardOp<T> operator+ (const BTypeName<T>& a, const BackwardOp<T>& b)
{
	ASSERT(a.m_val && b.m_val);
	return ADDBackwardOp<T>(a,b);
}

template <class T>
INLINE1 ADDBackwardOp<T> operator+ (const BTypeName<T>& a, const BTypeName<T>& b)
{
	ASSERT(a.m_val && b.m_val);
	return ADDBackwardOp<T>(a,b);
}

template <class T>
INLINE1 ADDBackwardOp3<T,typename Op<T>::Base> operator+ (const typename Op<T>::Base& a, const BackwardOp<T>& b)
{
	ASSERT(b.m_val);
	return ADDBackwardOp3<T,typename Op<T>::Base>(a,b);
}

template <class T>
INLINE1 ADDBackwardOp3<T,typename Op<T>::Base> operator+ (const typename Op<T>::Base& a, const BTypeName<T>& b)
{
	ASSERT(b.m_val);
	return ADDBackwardOp3<T,typename Op<T>::Base>(a,b);
}

template <class T>
INLINE1 ADDBackwardOp4<T,typename Op<T>::Base> operator+ (const BackwardOp<T>& a, const typename Op<T>::Base& b)
{
	ASSERT(a.m_val);
	return ADDBackwardOp4<T,typename Op<T>::Base>(a,b);
}

template <class T>
INLINE1 ADDBackwardOp4<T,typename Op<T>::Base> operator+ (const BTypeName<T>& a, const typename Op<T>::Base& b)
{
	ASSERT(a.m_val);
	return ADDBackwardOp4<T,typename Op<T>::Base>(a,b);
}
/*
template <class T>
INLINE1 ADDBackwardOp3<T,T> operator+ (const T& a, const BackwardOp<T>& b)
{
	ASSERT(b.m_val);
	return ADDBackwardOp3<T,T>(a,b);
}

template <class T>
INLINE1 ADDBackwardOp3<T,T> operator+ (const T& a, const BTypeName<T>& b)
{
	ASSERT(b.m_val);
	return ADDBackwardOp3<T,T>(a,b);
}

template <class T>
INLINE1 ADDBackwardOp4<T,T> operator+ (const BackwardOp<T>& a, const T& b)
{
	ASSERT(a.m_val);
	return ADDBackwardOp4<T,T>(a,b);
}

template <class T>
INLINE1 ADDBackwardOp4<T,T> operator+ (const BTypeName<T>& a, const T& b)
{
	ASSERT(a.m_val);
	return ADDBackwardOp4<T,T>(a,b);
}
*/
template <class T>
INLINE1 SUBBackwardOp<T> operator- (const BackwardOp<T>& a, const BackwardOp<T>& b)
{
	ASSERT(a.m_val && b.m_val);
	return SUBBackwardOp<T>(a,b);
}

template <class T>
INLINE1 SUBBackwardOp<T> operator- (const BackwardOp<T>& a, const BTypeName<T>& b)
{
	ASSERT(a.m_val && b.m_val);
	return SUBBackwardOp<T>(a,b);
}

template <class T>
INLINE1 SUBBackwardOp<T> operator- (const BTypeName<T>& a, const BackwardOp<T>& b)
{
	ASSERT(a.m_val && b.m_val);
	return SUBBackwardOp<T>(a,b);
}

template <class T>
INLINE1 SUBBackwardOp<T> operator- (const BTypeName<T>& a, const BTypeName<T>& b)
{
	ASSERT(a.m_val && b.m_val);
	return SUBBackwardOp<T>(a,b);
}

template <class T>
INLINE1 SUBBackwardOp3<T,typename Op<T>::Base> operator- (const typename Op<T>::Base& a, const BackwardOp<T>& b)
{
	ASSERT(b.m_val);
	return SUBBackwardOp3<T,typename Op<T>::Base>(a,b);
}

template <class T>
INLINE1 SUBBackwardOp3<T,typename Op<T>::Base> operator- (const typename Op<T>::Base& a, const BTypeName<T>& b)
{
	ASSERT(b.m_val);
	return SUBBackwardOp3<T,typename Op<T>::Base>(a,b);
}

template <class T>
INLINE1 SUBBackwardOp4<T,typename Op<T>::Base> operator- (const BackwardOp<T>& a, const typename Op<T>::Base& b)
{
	ASSERT(a.m_val);
	return SUBBackwardOp4<T,typename Op<T>::Base>(a,b);
}

template <class T>
INLINE1 SUBBackwardOp4<T,typename Op<T>::Base> operator- (const BTypeName<T>& a, const typename Op<T>::Base& b)
{
	ASSERT(a.m_val);
	return SUBBackwardOp4<T,typename Op<T>::Base>(a,b);
}
/*
template <class T>
INLINE1 SUBBackwardOp3<T,T> operator- (const T& a, const BackwardOp<T>& b)
{
	ASSERT(b.m_val);
	return SUBBackwardOp3<T,T>(a,b);
}

template <class T>
INLINE1 SUBBackwardOp3<T,T> operator- (const T& a, const BTypeName<T>& b)
{
	ASSERT(b.m_val);
	return SUBBackwardOp3<T,T>(a,b);
}

template <class T>
INLINE1 SUBBackwardOp4<T,T> operator- (const BackwardOp<T>& a, const T& b)
{
	ASSERT(a.m_val);
	return SUBBackwardOp4<T,T>(a,b);
}

template <class T>
INLINE1 SUBBackwardOp4<T,T> operator- (const BTypeName<T>& a, const T& b)
{
	ASSERT(a.m_val);
	return SUBBackwardOp4<T,T>(a,b);
}
*/
template <class T>
INLINE1 MULBackwardOp<T> operator* (const BackwardOp<T>& a, const BackwardOp<T>& b)
{
	ASSERT(a.m_val && b.m_val);
	return MULBackwardOp<T>(a,b);
}

template <class T>
INLINE1 MULBackwardOp<T> operator* (const BackwardOp<T>& a, const BTypeName<T>& b)
{
	return MULBackwardOp<T>(a,b);
}

template <class T>
INLINE1 MULBackwardOp<T> operator* (const BTypeName<T>& a, const BackwardOp<T>& b)
{
	ASSERT(a.m_val && b.m_val);
	return MULBackwardOp<T>(a,b);
}

template <class T>
INLINE1 MULBackwardOp<T> operator* (const BTypeName<T>& a, const BTypeName<T>& b)
{
	ASSERT(a.m_val && b.m_val);
	return MULBackwardOp<T>(a,b);
}

template <class T>
INLINE1 MULBackwardOp3<T,typename Op<T>::Base> operator* (const typename Op<T>::Base& a, const BackwardOp<T>& b)
{
	ASSERT(b.m_val);
	return MULBackwardOp3<T,typename Op<T>::Base>(a,b);
}

template <class T>
INLINE1 MULBackwardOp3<T,typename Op<T>::Base> operator* (const typename Op<T>::Base& a, const BTypeName<T>& b)
{
	ASSERT(b.m_val);
	return MULBackwardOp3<T,typename Op<T>::Base>(a,b);
}

template <class T>
INLINE1 MULBackwardOp4<T,typename Op<T>::Base> operator* (const BackwardOp<T>& a, const typename Op<T>::Base& b)
{
	ASSERT(a.m_val);
	return MULBackwardOp4<T,typename Op<T>::Base>(a,b);
}

template <class T>
INLINE1 MULBackwardOp4<T,typename Op<T>::Base> operator* (const BTypeName<T>& a, const typename Op<T>::Base& b)
{
	ASSERT(a.m_val);
	return MULBackwardOp4<T,typename Op<T>::Base>(a,b);
}
/*
template <class T>
INLINE1 MULBackwardOp3<T,T> operator* (const T& a, const BackwardOp<T>& b)
{
	ASSERT(b.m_val);
	return MULBackwardOp3<T,T>(a,b);
}

template <class T>
INLINE1 MULBackwardOp3<T,T> operator* (const T& a, const BTypeName<T>& b)
{
	ASSERT(b.m_val);
	return MULBackwardOp3<T,T>(a,b);
}

template <class T>
INLINE1 MULBackwardOp4<T,T> operator* (const BackwardOp<T>& a, const T& b)
{
	ASSERT(a.m_val);
	return MULBackwardOp4<T,T>(a,b);
}

template <class T>
INLINE1 MULBackwardOp4<T,T> operator* (const BTypeName<T>& a, const T& b)
{
	ASSERT(a.m_val);
	return MULBackwardOp4<T,T>(a,b);
}
*/
template <class T>
INLINE1 DIVBackwardOp<T> operator/ (const BackwardOp<T>& a, const BackwardOp<T>& b)
{
	ASSERT(a.m_val && b.m_val);
	return DIVBackwardOp<T>(a,b);
}

template <class T>
INLINE1 DIVBackwardOp<T> operator/ (const BackwardOp<T>& a, const BTypeName<T>& b)
{
	ASSERT(a.m_val && b.m_val);
	return DIVBackwardOp<T>(a,b);
}

template <class T>
INLINE1 DIVBackwardOp<T> operator/ (const BTypeName<T>& a, const BackwardOp<T>& b)
{
	ASSERT(a.m_val && b.m_val);
	return DIVBackwardOp<T>(a,b);
}

template <class T>
INLINE1 DIVBackwardOp<T> operator/ (const BTypeName<T>& a, const BTypeName<T>& b)
{
	ASSERT(a.m_val && b.m_val);
	return DIVBackwardOp<T>(a,b);
}

template <class T>
INLINE1 DIVBackwardOp3<T,typename Op<T>::Base> operator/ (const typename Op<T>::Base& a, const BackwardOp<T>& b)
{
	ASSERT(b.m_val);
	return DIVBackwardOp3<T,typename Op<T>::Base>(a,b);
}

template <class T>
INLINE1 DIVBackwardOp3<T,typename Op<T>::Base> operator/ (const typename Op<T>::Base& a, const BTypeName<T>& b)
{
	ASSERT(b.m_val);
	return DIVBackwardOp3<T,typename Op<T>::Base>(a,b);
}

template <class T>
INLINE1 DIVBackwardOp4<T,typename Op<T>::Base> operator/ (const BackwardOp<T>& a, const typename Op<T>::Base& b)
{
	ASSERT(a.m_val);
	return DIVBackwardOp4<T,typename Op<T>::Base>(a,b);
}

template <class T>
INLINE1 DIVBackwardOp4<T,typename Op<T>::Base> operator/ (const BTypeName<T>& a, const typename Op<T>::Base& b)
{
	ASSERT(a.m_val);
	return DIVBackwardOp4<T,typename Op<T>::Base>(a,b);
}
/*
template <class T>
INLINE1 DIVBackwardOp3<T,T> operator/ (const T& a, const BackwardOp<T>& b)
{
	ASSERT(b.m_val);
	return DIVBackwardOp3<T,typename Op<T>::Base>(a,b);
}

template <class T>
INLINE1 DIVBackwardOp3<T,T> operator/ (const T& a, const BTypeName<T>& b)
{
	ASSERT(b.m_val);
	return DIVBackwardOp3<T,T>(a,b);
}

template <class T>
INLINE1 DIVBackwardOp4<T,T> operator/ (const BackwardOp<T>& a, const T& b)
{
	ASSERT(a.m_val);
	return DIVBackwardOp4<T,T>(a,b);
}

template <class T>
INLINE1 DIVBackwardOp4<T,T> operator/ (const BTypeName<T>& a, const T& b)
{
	ASSERT(a.m_val);
	return DIVBackwardOp4<T,T>(a,b);
}
*/
template <class T>
INLINE1 POWBackwardOp<T> pow (const BackwardOp<T>& a, const BackwardOp<T>& b)
{
	ASSERT(a.m_val && b.m_val);
	return POWBackwardOp<T>(a,b);
}

template <class T>
INLINE1 POWBackwardOp<T> pow (const BackwardOp<T>& a, const BTypeName<T>& b)
{
	ASSERT(a.m_val && b.m_val);
	return POWBackwardOp<T>(a,b);
}

template <class T>
INLINE1 POWBackwardOp<T> pow (const BTypeName<T>& a, const BackwardOp<T>& b)
{
	ASSERT(a.m_val && b.m_val);
	return POWBackwardOp<T>(a,b);
}

template <class T>
INLINE1 POWBackwardOp<T> pow (const BTypeName<T>& a, const BTypeName<T>& b)
{
	ASSERT(a.m_val && b.m_val);
	return POWBackwardOp<T>(a,b);
}

template <class T>
INLINE1 POWBackwardOp3<T> pow (const typename Op<T>::Base& a, const BackwardOp<T>& b)
{
	ASSERT(b.m_val);
	return POWBackwardOp3<T>(a,b);
}

template <class T>
INLINE1 POWBackwardOp3<T> pow (const typename Op<T>::Base& a, const BTypeName<T>& b)
{
	ASSERT(b.m_val);
	return POWBackwardOp3<T>(a,b);
}

template <class T>
INLINE1 POWBackwardOp4<T> pow (const BackwardOp<T>& a, const typename Op<T>::Base& b)
{
	ASSERT(a.m_val);
	return POWBackwardOp4<T>(a,b);
}

template <class T>
INLINE1 POWBackwardOp4<T> pow (const BTypeName<T>& a, const typename Op<T>::Base& b)
{
	ASSERT(a.m_val);
	return POWBackwardOp4<T>(a,b);
}

template <class T>
INLINE1 UADDBackwardOp<T> operator+ (const BackwardOp<T>& a)
{
	ASSERT(a.m_val);
	return UADDBackwardOp<T>(a);
}

template <class T>
INLINE1 USUBBackwardOp<T> operator- (const BackwardOp<T>& a)
{
	ASSERT(a.m_val);
	return USUBBackwardOp<T>(a);
}

template <class T>
INLINE1 POWopN<T> pow (const BackwardOp<T>& a, int ib)
{
	ASSERT(a.m_val);
	return POWopN<T>(a,ib);
}

template <class T>
INLINE1 POWopN<T> pow (const BTypeName<T>& a, int ib)
{
	ASSERT(a.m_val);
	return POWopN<T>(a,ib);
}

template <class T>
INLINE1 SQRBackwardOp<T> sqr (const BackwardOp<T>& a)
{
	ASSERT(a.m_val);
	return SQRBackwardOp<T>(a);
}

template <class T>
INLINE1 EXPBackwardOp<T> exp (const BackwardOp<T>& a)
{
	ASSERT(a.m_val);
	return EXPBackwardOp<T>(a);
}

template <class T>
INLINE1 LOGBackwardOp<T> log (const BackwardOp<T>& a)
{
	ASSERT(a.m_val);
	return LOGBackwardOp<T>(a);
}

template <class T>
INLINE1 SQRTBackwardOp<T> sqrt (const BackwardOp<T>& a)
{
	ASSERT(a.m_val);
	return SQRTBackwardOp<T>(a);
}

template <class T>
INLINE1 SINBackwardOp<T> sin (const BackwardOp<T>& a)
{
	ASSERT(a.m_val);
	return SINBackwardOp<T>(a);
}

template <class T>
INLINE1 COSBackwardOp<T> cos (const BackwardOp<T>& a)
{
	ASSERT(a.m_val);
	return COSBackwardOp<T>(a);
}

template <class T>
INLINE1 TANBackwardOp<T> tan (const BackwardOp<T>& a)
{
	ASSERT(a.m_val);
	return TANBackwardOp<T>(a);
}

template <class T>
INLINE1 ASINBackwardOp<T> asin (const BackwardOp<T>& a)
{
	ASSERT(a.m_val);
	return ASINBackwardOp<T>(a);
}

template <class T>
INLINE1 ACOSBackwardOp<T> acos (const BackwardOp<T>& a)
{
	ASSERT(a.m_val);
	return ACOSBackwardOp<T>(a);
}

template <class T>
INLINE1 ATANBackwardOp<T> atan (const BackwardOp<T>& a)
{
	ASSERT(a.m_val);
	return ATANBackwardOp<T>(a);
}

template <class T> template <class U>
INLINE1 BTypeName<T>& BTypeName<T>::operator += (const U& x)
{
	return *this = (*this)+x;
}

template <class T> template <class U>
INLINE1 BTypeName<T>& BTypeName<T>::operator -= (const U& x)
{
	return *this = (*this)-x;
}

template <class T> template <class U>
INLINE1 BTypeName<T>& BTypeName<T>::operator *= (const U& x)
{
	return *this = (*this)*x;
}

template <class T> template <class U>
INLINE1 BTypeName<T>& BTypeName<T>::operator /= (const U& x)
{
	return *this = (*this)/x;
}

template <typename U> struct Op< BTypeName<U> >
{
	typedef BTypeName<U> T;
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
