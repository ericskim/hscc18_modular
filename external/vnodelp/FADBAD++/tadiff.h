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
// COPYRIGHT NOTICE
// ***************************************************************

#ifndef _TADIFF_H
#define _TADIFF_H

#include <algorithm>

#ifndef MaxLength
#define MaxLength 40
#endif

#include "fadbad.h"

template <typename U, int N>
class TValues
{
	unsigned int m_n;
	U m_val[N];
public:	
	TValues():m_n(0){std::fill(m_val,m_val+N,Op<U>::myZero());}
	template <typename V> explicit TValues(const V& val):m_n(1){m_val[0]=val;std::fill(m_val+1,m_val+N,Op<U>::myZero());}
	U& operator[](const unsigned int i) { return m_val[i]; }
	const U& operator[](const unsigned int i) const { return m_val[i]; }
	const unsigned int length() const { return m_n; }
	unsigned int& length() { return m_n; }
	void reset(){ m_n=0; }
};

template <typename U, int N>
class TTypeNameHV // Heap Value
{
	TValues<U,N> m_val;
	mutable unsigned int m_rc;
protected:
	virtual ~TTypeNameHV(){}
public:
	TTypeNameHV():m_rc(0){}
	template <typename V> explicit TTypeNameHV(const V& val):m_val(val),m_rc(0){}
	const U& val(const unsigned int i) const { return m_val[i]; }
	U& val(const unsigned int i) { return m_val[i]; }
	const unsigned int length() const { return m_val.length(); }
	unsigned int& length() { return m_val.length(); }
	void decRef(TTypeNameHV<U,N>*& pTTypeNameHV) const { if (--m_rc==0) { delete this; pTTypeNameHV=0;} }
	void incRef() const {++m_rc;}

	virtual void reset(){m_val.reset();}
	virtual unsigned int eval(const unsigned int k){return k+1;}
};

template <typename U, int N=MaxLength>
class TTypeName
{
private:
	struct SV // Stack Value refers to reference-counted Heap Value:
	{
		mutable TTypeNameHV<U,N>* m_pTTypeNameHV;
		SV(TTypeNameHV<U,N>* pTTypeNameHV):m_pTTypeNameHV(pTTypeNameHV){ m_pTTypeNameHV->incRef(); }
		SV(const typename TTypeName<U,N>::SV& sv):m_pTTypeNameHV(sv.m_pTTypeNameHV){ m_pTTypeNameHV->incRef(); }
		~SV(){ m_pTTypeNameHV->decRef(m_pTTypeNameHV); }
		TTypeNameHV<U,N>* getTTypeNameHV() const { return m_pTTypeNameHV; }
		void setTTypeNameHV(TTypeNameHV<U,N>* pTTypeNameHV) 
		{ 
			if (m_pTTypeNameHV!=pTTypeNameHV) 
			{ 
				m_pTTypeNameHV->decRef(m_pTTypeNameHV);
				m_pTTypeNameHV=pTTypeNameHV;
				m_pTTypeNameHV->incRef();
			}
		}
		const U& val() const { return m_pTTypeNameHV->val(0); }
		const U& val(const unsigned int i) const { return m_pTTypeNameHV->val(i); }
		const unsigned int length() const { return m_pTTypeNameHV->length(); }
		unsigned int& length() { return m_pTTypeNameHV->length(); }	
		U& val(const unsigned int i) { return m_pTTypeNameHV->val(i); }

		void reset(){m_pTTypeNameHV->reset();}
		unsigned int eval(const unsigned int i){return m_pTTypeNameHV->eval(i);}
	} m_sv;
public:
	TTypeName():m_sv(new TTypeNameHV<U,N>()){}
	TTypeName(TTypeNameHV<U,N>* pTTypeNameHV):m_sv(pTTypeNameHV){}
	explicit TTypeName(const typename TTypeName<U,N>::SV& sv):m_sv(sv){}
	template <typename V> /*explicit*/ TTypeName(const V& val):m_sv(new TTypeNameHV<U,N>(val)){m_sv.length()=N;}
	TTypeName<U,N>& operator=(const TTypeName<U,N>& val) { m_sv.setTTypeNameHV(val.m_sv.getTTypeNameHV()); return *this; }
	template <typename V> TTypeName<U,N>& operator=(const V& val) 
	{ 
		m_sv.setTTypeNameHV(new TTypeNameHV<U,N>(val));
		m_sv.length()=N;
		return *this; 
	}
	TTypeNameHV<U,N>* getTTypeNameHV() const { return m_sv.getTTypeNameHV(); }
	void setTTypeNameHV(const TTypeNameHV<U,N>* pTTypeNameHV) { m_sv.setTTypeNameHV(pTTypeNameHV); }
	const U& val() const { return m_sv.val(); }
	const unsigned int length() const { return m_sv.length(); }	
	const U& operator[](const unsigned int i) const { return m_sv.val(i); }
	U& operator[](const unsigned int i) { if (i>=m_sv.length()) m_sv.length()=i+1; return m_sv.val(i);}
	
	TTypeName<U,N>& operator+=(const TTypeName<U,N>& val);
	TTypeName<U,N>& operator-=(const TTypeName<U,N>& val);
	TTypeName<U,N>& operator*=(const TTypeName<U,N>& val);
	TTypeName<U,N>& operator/=(const TTypeName<U,N>& val);
	template <typename V> TTypeName<U,N>& operator+=(const V& val);
	template <typename V> TTypeName<U,N>& operator-=(const V& val);
	template <typename V> TTypeName<U,N>& operator*=(const V& val);
	template <typename V> TTypeName<U,N>& operator/=(const V& val);	

	void reset(){m_sv.reset();}
	unsigned int eval(const unsigned int i){return m_sv.eval(i);}
};

template <typename U, int N> bool operator==(const TTypeName<U,N>& val1, const TTypeName<U,N>& val2) { return val1.val()==val2.val(); }
template <typename U, int N> bool operator!=(const TTypeName<U,N>& val1, const TTypeName<U,N>& val2) { return val1.val()!=val2.val(); }
template <typename U, int N> bool operator<(const TTypeName<U,N>& val1, const TTypeName<U,N>& val2) { return val1.val()<val2.val(); }
template <typename U, int N> bool operator<=(const TTypeName<U,N>& val1, const TTypeName<U,N>& val2) { return val1.val()<=val2.val(); }
template <typename U, int N> bool operator>(const TTypeName<U,N>& val1, const TTypeName<U,N>& val2) { return val1.val()>val2.val(); }
template <typename U, int N> bool operator>=(const TTypeName<U,N>& val1, const TTypeName<U,N>& val2) { return val1.val()>=val2.val(); }
template <typename U, int N, typename V> bool operator==(const TTypeName<U,N>& val1, const V& val2) { return val1.val()==val2; }
template <typename U, int N, typename V> bool operator==(const V& val1, const TTypeName<U,N>& val2) { return val1==val2.val(); }
template <typename U, int N, typename V> bool operator!=(const TTypeName<U,N>& val1, const V& val2) { return val1.val()!=val2; }
template <typename U, int N, typename V> bool operator!=(const V& val1, const TTypeName<U,N>& val2) { return val1!=val2.val(); }
template <typename U, int N, typename V> bool operator<(const TTypeName<U,N>& val1, const V& val2) { return val1.val()<val2; }
template <typename U, int N, typename V> bool operator<(const V& val1, const TTypeName<U,N>& val2) { return val1<val2.val(); }
template <typename U, int N, typename V> bool operator<=(const TTypeName<U,N>& val1, const V& val2) { return val1.val()<=val2; }
template <typename U, int N, typename V> bool operator<=(const V& val1, const TTypeName<U,N>& val2) { return val1<=val2.val(); }
template <typename U, int N, typename V> bool operator>(const TTypeName<U,N>& val1, const V& val2) { return val1.val()>val2; }
template <typename U, int N, typename V> bool operator>(const V& val1, const TTypeName<U,N>& val2) { return val1>val2.val(); }
template <typename U, int N, typename V> bool operator>=(const TTypeName<U,N>& val1, const V& val2) { return val1.val()>=val2; }
template <typename U, int N, typename V> bool operator>=(const V& val1, const TTypeName<U,N>& val2) { return val1>=val2.val(); }

// Binary operator base class:

template <typename U, int N>
class BinTTypeNameHV : public TTypeNameHV<U,N>
{
	TTypeNameHV<U,N>* m_pOp1;
	TTypeNameHV<U,N>* m_pOp2;
public:
	BinTTypeNameHV(const U& val, TTypeNameHV<U,N>* pOp1, TTypeNameHV<U,N>* pOp2):TTypeNameHV<U,N>(val),m_pOp1(pOp1),m_pOp2(pOp2)
	{
		m_pOp1->incRef();m_pOp2->incRef();
	}
	BinTTypeNameHV(TTypeNameHV<U,N>* pOp1, TTypeNameHV<U,N>* pOp2):TTypeNameHV<U,N>(),m_pOp1(pOp1),m_pOp2(pOp2)
	{
		m_pOp1->incRef();m_pOp2->incRef();
	}
	virtual ~BinTTypeNameHV()
	{
		m_pOp1->decRef(m_pOp1);m_pOp2->decRef(m_pOp2);
	}
	TTypeNameHV<U,N>* op1() { return m_pOp1; }
	TTypeNameHV<U,N>* op2() { return m_pOp2; }

	unsigned int op1Eval(const unsigned int k){return this->op1()->eval(k);}
	unsigned int op2Eval(const unsigned int k){return this->op2()->eval(k);}
	const U& op1Val(const unsigned int k) {return this->op1()->val(k);}
	const U& op2Val(const unsigned int k) {return this->op2()->val(k);}
	void reset(){op1()->reset();op2()->reset();TTypeNameHV<U,N>::reset();}
};

// Unary operator base class:

template <typename U, int N>
class UnTTypeNameHV : public TTypeNameHV<U,N>
{
	TTypeNameHV<U,N>* m_pOp;
public:
	UnTTypeNameHV(const U& val, TTypeNameHV<U,N>* pOp):TTypeNameHV<U,N>(val),m_pOp(pOp)
	{
		m_pOp->incRef();
	}
	UnTTypeNameHV(TTypeNameHV<U,N>* pOp):TTypeNameHV<U,N>(),m_pOp(pOp)
	{
		m_pOp->incRef();
	}
	virtual ~UnTTypeNameHV()
	{
		m_pOp->decRef(m_pOp);
	}
	TTypeNameHV<U,N>* op() { return m_pOp; }

	unsigned int opEval(const unsigned int k){return this->op()->eval(k);}
	const U& opVal(const unsigned int k) {return this->op()->val(k);}
	void reset(){op()->reset();TTypeNameHV<U,N>::reset();}
};

// ADDITION:

template <typename U, int N>
struct ADD : public BinTTypeNameHV<U,N>
{
	ADD(const U& val, TTypeNameHV<U,N>* pOp1, TTypeNameHV<U,N>* pOp2):BinTTypeNameHV<U,N>(val,pOp1,pOp2){}
	ADD(TTypeNameHV<U,N>* pOp1, TTypeNameHV<U,N>* pOp2):BinTTypeNameHV<U,N>(pOp1,pOp2){}
	unsigned int eval(const unsigned int k)
	{
		unsigned int l=std::min(this->op1Eval(k),this->op2Eval(k));
		for(unsigned int i=this->length();i<l;++i) this->val(i)=this->op1Val(i)+this->op2Val(i);
		return this->length()=l;
	}
private:
	void operator=(const ADD<U,N>&){} // not allowed
};
template <typename U, int N, typename V>
struct ADD1 : public UnTTypeNameHV<U,N>
{
	const V m_a;
	ADD1(const U& val, const V& a, TTypeNameHV<U,N>* pOp2):UnTTypeNameHV<U,N>(val,pOp2),m_a(a){}
	ADD1(const V& a, TTypeNameHV<U,N>* pOp2):UnTTypeNameHV<U,N>(pOp2),m_a(a){}
	unsigned int eval(const unsigned int k)
	{
		unsigned int l=this->opEval(k);
		if (0==this->length()) { this->val(0)=m_a+this->opVal(0); this->length()=1; }
		for(unsigned int i=this->length();i<l;++i) this->val(i)=this->opVal(i);
		return this->length()=l;
	}
private:
	void operator=(const ADD1<U,N,V>&){} // not allowed
};
template <typename U, int N, typename V>
struct ADD2 : public UnTTypeNameHV<U,N>
{
	const V m_b;
	ADD2(const U& val, TTypeNameHV<U,N>* pOp1, const V& b):UnTTypeNameHV<U,N>(val,pOp1),m_b(b){}
	ADD2(TTypeNameHV<U,N>* pOp1, const V& b):UnTTypeNameHV<U,N>(pOp1),m_b(b){}
	unsigned int eval(const unsigned int k)
	{
		unsigned int l=this->opEval(k);
		if (0==this->length()) { this->val(0)=this->opVal(0)+m_b; this->length()=1; }
		for(unsigned int i=this->length();i<l;++i) this->val(i)=this->opVal(i);
		return this->length()=l;
	}
private:
	void operator=(const ADD2<U,N,V>&){} // not allowed
};
template <typename U, int N>
TTypeName<U,N> operator+(const TTypeName<U,N>& val1, const TTypeName<U,N>& val2)
{
	TTypeNameHV<U,N>* pHV=val1.length()>0 && val2.length()>0 ?
		new ADD<U,N>(val1.val()+val2.val(),val1.getTTypeNameHV(),val2.getTTypeNameHV()):
		new ADD<U,N>(val1.getTTypeNameHV(),val2.getTTypeNameHV());
	return TTypeName<U,N>(pHV);
}
template <typename U, int N, typename V>
TTypeName<U,N> operator+(const V& a, const TTypeName<U,N>& val2)
{
	TTypeNameHV<U,N>* pHV=val2.length()>0 ?
		new ADD1<U,N,V>(a+val2.val(), a, val2.getTTypeNameHV()):
		new ADD1<U,N,V>(a, val2.getTTypeNameHV());
	return TTypeName<U,N>(pHV);
}
template <typename U, int N, typename V>
TTypeName<U,N> operator+(const TTypeName<U,N>& val1, const V& b)
{
	TTypeNameHV<U,N>* pHV=val1.length()>0?
		new ADD2<U,N,V>(val1.val()+b, val1.getTTypeNameHV(), b):
		new ADD2<U,N,V>(val1.getTTypeNameHV(), b);
	return TTypeName<U,N>(pHV);
}

// SUBTRACTION:

template <typename U, int N>
struct SUB : public BinTTypeNameHV<U,N>
{
	SUB(const U& val, TTypeNameHV<U,N>* pOp1, TTypeNameHV<U,N>* pOp2):BinTTypeNameHV<U,N>(val,pOp1,pOp2){}
	SUB(TTypeNameHV<U,N>* pOp1, TTypeNameHV<U,N>* pOp2):BinTTypeNameHV<U,N>(pOp1,pOp2){}
	unsigned int eval(const unsigned int k)
	{
		unsigned int l=std::min(this->op1Eval(k),this->op2Eval(k));
		for(unsigned int i=this->length();i<l;++i) this->val(i)=this->op1Val(i)-this->op2Val(i);
		return this->length()=l;
	}
private:
	void operator=(const SUB<U,N>&){} // not allowed
};
template <typename U, int N, typename V>
struct SUB1 : public UnTTypeNameHV<U,N>
{
	const V m_a;
	SUB1(const U& val, const V& a, TTypeNameHV<U,N>* pOp2):UnTTypeNameHV<U,N>(val,pOp2),m_a(a){}
	SUB1(const V& a, TTypeNameHV<U,N>* pOp2):UnTTypeNameHV<U,N>(pOp2),m_a(a){}
	unsigned int eval(const unsigned int k)
	{
		unsigned int l=this->opEval(k);
		if (0==this->length()) { this->val(0)=m_a-this->opVal(0); this->length()=1; }
		for(unsigned int i=this->length();i<l;++i) this->val(i)=-this->opVal(i);
		return this->length()=l;
	}
private:
	void operator=(const SUB1<U,N,V>&){} // not allowed
};
template <typename U, int N, typename V>
struct SUB2 : public UnTTypeNameHV<U,N>
{
	const V m_b;
	SUB2(const U& val, TTypeNameHV<U,N>* pOp1, const V& b):UnTTypeNameHV<U,N>(val,pOp1),m_b(b){}
	SUB2(TTypeNameHV<U,N>* pOp1, const V& b):UnTTypeNameHV<U,N>(pOp1),m_b(b){}
	unsigned int eval(const unsigned int k)
	{
		unsigned int l=this->opEval(k);
		if (0==this->length()) { this->val(0)=this->opVal(0)-m_b; this->length()=1; }
		for(unsigned int i=this->length();i<l;++i) this->val(i)=this->opVal(i);
		return this->length()=l;
	}
private:
	void operator=(const SUB2<U,N,V>&){} // not allowed
};
template <typename U, int N>
TTypeName<U,N> operator-(const TTypeName<U,N>& val1, const TTypeName<U,N>& val2)
{ 
	TTypeNameHV<U,N>* pHV=val1.length()>0 && val2.length()>0 ?
		new SUB<U,N>(val1.val()-val2.val(),val1.getTTypeNameHV(),val2.getTTypeNameHV()):
		new SUB<U,N>(val1.getTTypeNameHV(),val2.getTTypeNameHV());
	return TTypeName<U,N>(pHV);
}
template <typename U, int N, typename V>
TTypeName<U,N> operator-(const V& a, const TTypeName<U,N>& val2)
{
	TTypeNameHV<U,N>* pHV=val2.length()>0 ?
		new SUB1<U,N,V>(a-val2.val(), a, val2.getTTypeNameHV()):
		new SUB1<U,N,V>(a, val2.getTTypeNameHV());
	return TTypeName<U,N>(pHV);
}
template <typename U, int N, typename V>
TTypeName<U,N> operator-(const TTypeName<U,N>& val1, const V& b)
{
	TTypeNameHV<U,N>* pHV=val1.length()>0 ?
		new SUB2<U,N,V>(val1.val()-b, val1.getTTypeNameHV(), b):
		new SUB2<U,N,V>(val1.getTTypeNameHV(), b);
	return TTypeName<U,N>(pHV);
}

// MULTIPLICATION:

template <typename U, int N>
struct MUL : public BinTTypeNameHV<U,N>
{
	MUL(const U& val, TTypeNameHV<U,N>* pOp1, TTypeNameHV<U,N>* pOp2):BinTTypeNameHV<U,N>(val,pOp1,pOp2){}
	MUL(TTypeNameHV<U,N>* pOp1, TTypeNameHV<U,N>* pOp2):BinTTypeNameHV<U,N>(pOp1,pOp2){}
	unsigned int eval(const unsigned int k)
	{
		unsigned int l=std::min(this->op1Eval(k),this->op2Eval(k));
		for(unsigned int i=this->length();i<l;++i)
		{
			this->val(i)=Op<U>::myZero();
			for(unsigned int j=0;j<=i;++j)
				Op<U>::myCadd(this->val(i),this->op1Val(j)*this->op2Val(i-j));
		}
		return this->length()=l;
	}
private:
	void operator=(const MUL<U,N>&){} // not allowed
};
template <typename U, int N, typename V>
struct MUL1 : public UnTTypeNameHV<U,N>
{
	const V m_a;
	MUL1(const U& val, const V& a, TTypeNameHV<U,N>* pOp2):UnTTypeNameHV<U,N>(val,pOp2),m_a(a){}
	MUL1(const V& a, TTypeNameHV<U,N>* pOp2):UnTTypeNameHV<U,N>(pOp2),m_a(a){}
	unsigned int eval(const unsigned int k)
	{
		unsigned int l=this->opEval(k);
		for(unsigned int i=this->length();i<l;++i) this->val(i)=m_a*this->opVal(i);
		return this->length()=l;
	}
private:
	void operator=(const MUL1<U,N,V>&){} // not allowed
};
template <typename U, int N, typename V>
struct MUL2 : public UnTTypeNameHV<U,N>
{
	const V m_b;
	MUL2(const U& val, TTypeNameHV<U,N>* pOp1, const V& b):UnTTypeNameHV<U,N>(val,pOp1),m_b(b){}
	MUL2(TTypeNameHV<U,N>* pOp1, const V& b):UnTTypeNameHV<U,N>(pOp1),m_b(b){}
	unsigned int eval(const unsigned int k)
	{
		unsigned int l=this->opEval(k);
		for(unsigned int i=this->length();i<l;++i) this->val(i)=this->opVal(i)*m_b;
		return this->length()=l;
	}
private:
	void operator=(const MUL2<U,N,V>&){} // not allowed
};
template <typename U, int N>
TTypeName<U,N> operator*(const TTypeName<U,N>& val1, const TTypeName<U,N>& val2)
{
	TTypeNameHV<U,N>* pHV=val1.length()>0 && val2.length()>0 ?
		new MUL<U,N>(val1.val()*val2.val(),val1.getTTypeNameHV(),val2.getTTypeNameHV()):
		new MUL<U,N>(val1.getTTypeNameHV(),val2.getTTypeNameHV());
	return TTypeName<U,N>(pHV);
}
template <typename U, int N, typename V>
TTypeName<U,N> operator*(const V& a, const TTypeName<U,N>& val2)
{
	TTypeNameHV<U,N>* pHV=val2.length()>0 ?
		new MUL1<U,N,V>(a*val2.val(), a, val2.getTTypeNameHV()):
		new MUL1<U,N,V>(a, val2.getTTypeNameHV());
	return TTypeName<U,N>(pHV);
}
template <typename U, int N, typename V>
TTypeName<U,N> operator*(const TTypeName<U,N>& val1, const V& b)
{
	TTypeNameHV<U,N>* pHV=val1.length()>0 ?
		new MUL2<U,N,V>(val1.val()*b, val1.getTTypeNameHV(), b):
		new MUL2<U,N,V>(val1.getTTypeNameHV(), b);
	return TTypeName<U,N>(pHV);
}

// DIVISION:

template <typename U, int N>
struct DIV : public BinTTypeNameHV<U,N>
{
	DIV(const U& val, TTypeNameHV<U,N>* pOp1, TTypeNameHV<U,N>* pOp2):BinTTypeNameHV<U,N>(val,pOp1,pOp2){}
	DIV(TTypeNameHV<U,N>* pOp1, TTypeNameHV<U,N>* pOp2):BinTTypeNameHV<U,N>(pOp1,pOp2){}
	unsigned int eval(const unsigned int k)
	{
		unsigned int l=std::min(this->op1Eval(k),this->op2Eval(k));
		for(unsigned int i=this->length();i<l;++i)
		{
			this->val(i)=this->op1Val(i);
			for(unsigned int j=1;j<=i;++j) Op<U>::myCsub(this->val(i),this->op2Val(j)*this->val(i-j));
			Op<U>::myCdiv(this->val(i), this->op2Val(0));
		}
		return this->length()=l;
	}
private:
	void operator=(const DIV<U,N>&){} // not allowed
};
template <typename U, int N, typename V>
struct DIV1 : public UnTTypeNameHV<U,N>
{
	const V m_a;
	DIV1(const U& val, const V& a, TTypeNameHV<U,N>* pOp2):UnTTypeNameHV<U,N>(val,pOp2),m_a(a){}
	DIV1(const V& a, TTypeNameHV<U,N>* pOp2):UnTTypeNameHV<U,N>(pOp2),m_a(a){}
	unsigned int eval(const unsigned int k)
	{
		unsigned int l=this->opEval(k);
		if (0==this->length()) { this->val(0)=m_a/this->opVal(0); this->length()=1; }
		for(unsigned int i=this->length();i<l;++i)
		{
			this->val(i)=Op<U>::myZero();
			for(unsigned int j=1;j<=i;++j) Op<U>::myCsub(this->val(i),this->opVal(j)*this->val(i-j));
			Op<U>::myCdiv(this->val(i), this->opVal(0));
		}
		return this->length()=l;
	}
private:
	void operator=(const DIV1<U,N,V>&){} // not allowed
};
template <typename U, int N, typename V>
struct DIV2 : public UnTTypeNameHV<U,N>
{
	const V m_b;
	DIV2(const U& val, TTypeNameHV<U,N>* pOp1, const V& b):UnTTypeNameHV<U,N>(val,pOp1),m_b(b){}
	DIV2(TTypeNameHV<U,N>* pOp1, const V& b):UnTTypeNameHV<U,N>(pOp1),m_b(b){}
	unsigned int eval(const unsigned int k)
	{
		unsigned int l=this->opEval(k);
		for(unsigned int i=this->length();i<l;++i) this->val(i)=this->opVal(i)/m_b;
		return this->length()=l;
	}
private:
	void operator=(const DIV2<U,N,V>&){} // not allowed
};
template <typename U, int N>
TTypeName<U,N> operator/(const TTypeName<U,N>& val1, const TTypeName<U,N>& val2)
{
	TTypeNameHV<U,N>* pHV=val1.length()>0 && val2.length()>0 ?
		new DIV<U,N>(val1.val()/val2.val(),val1.getTTypeNameHV(),val2.getTTypeNameHV()):
		new DIV<U,N>(val1.getTTypeNameHV(),val2.getTTypeNameHV());
	return TTypeName<U,N>(pHV);
}
template <typename U, int N, typename V>
TTypeName<U,N> operator/(const V& a, const TTypeName<U,N>& val2)
{
	TTypeNameHV<U,N>* pHV=val2.length()>0 ?
		new DIV1<U,N,V>(a/val2.val(), a, val2.getTTypeNameHV()):
		new DIV1<U,N,V>(a, val2.getTTypeNameHV());
	return TTypeName<U,N>(pHV);
}
template <typename U, int N, typename V>
TTypeName<U,N> operator/(const TTypeName<U,N>& val1, const V& b)
{
	TTypeNameHV<U,N>* pHV=val1.length()>0 ?
		new DIV2<U,N,V>(val1.val()/b, val1.getTTypeNameHV(), b):
		new DIV2<U,N,V>(val1.getTTypeNameHV(), b);
	return TTypeName<U,N>(pHV);
}

// COMPOUND ASSIGNMENTS:

template <typename U, int N> TTypeName<U,N>& TTypeName<U,N>::operator+=(const TTypeName<U,N>& val) { return (*this)=(*this)+val; }
template <typename U, int N> TTypeName<U,N>& TTypeName<U,N>::operator-=(const TTypeName<U,N>& val) { return (*this)=(*this)-val; } 
template <typename U, int N> TTypeName<U,N>& TTypeName<U,N>::operator*=(const TTypeName<U,N>& val) { return (*this)=(*this)*val; }
template <typename U, int N> TTypeName<U,N>& TTypeName<U,N>::operator/=(const TTypeName<U,N>& val) { return (*this)=(*this)/val; }
template <typename U, int N> template <typename V> TTypeName<U,N>& TTypeName<U,N>::operator+=(const V& val) { return (*this)=(*this)+val; }
template <typename U, int N> template <typename V> TTypeName<U,N>& TTypeName<U,N>::operator-=(const V& val) { return (*this)=(*this)-val; }
template <typename U, int N> template <typename V> TTypeName<U,N>& TTypeName<U,N>::operator*=(const V& val) { return (*this)=(*this)*val; }
template <typename U, int N> template <typename V> TTypeName<U,N>& TTypeName<U,N>::operator/=(const V& val) { return (*this)=(*this)/val; }

// UNARY MINUS

template <typename U, int N>
struct UMINUS : public UnTTypeNameHV<U,N>
{
	UMINUS(const U& val, TTypeNameHV<U,N>* pOp):UnTTypeNameHV<U,N>(val,pOp){}
	UMINUS(TTypeNameHV<U,N>* pOp):UnTTypeNameHV<U,N>(pOp){}
	unsigned int eval(const unsigned int k)
	{
		unsigned int l=this->opEval(k);
		for(unsigned int i=this->length();i<l;++i) this->val(i)=-this->opVal(i);
		return this->length()=l;
	}
private:
	void operator=(const UMINUS<U,N>&){} // not allowed
};

template <typename U, int N>
TTypeName<U,N> operator-(const TTypeName<U,N>& val)
{
	TTypeNameHV<U,N>* pHV=val.length()>0 ? 
		new UMINUS<U,N>(-val.val(),val.getTTypeNameHV()) : 
		new UMINUS<U,N>(val.getTTypeNameHV());
	return TTypeName<U,N>(pHV);
}

// UNARY PLUS

template <typename U, int N>
struct UPLUS : public UnTTypeNameHV<U,N>
{
	UPLUS(const U& val, TTypeNameHV<U,N>* pOp):UnTTypeNameHV<U,N>(val,pOp){}
	UPLUS(TTypeNameHV<U,N>* pOp):UnTTypeNameHV<U,N>(pOp){}
	unsigned int eval(const unsigned int k)
	{
		unsigned int l=this->opEval(k);
		for(unsigned int i=this->length();i<l;++i) this->val(i)=+this->opVal(i);
		return this->length()=l;
	}
private:
	void operator=(const UPLUS<U,N>&){} // not allowed
};

template <typename U, int N>
TTypeName<U,N> operator+(const TTypeName<U,N>& val)
{
	TTypeNameHV<U,N>* pHV=val.length()>0 ? 
		new UPLUS<U,N>(+val.val(),val.getTTypeNameHV()) : 
		new UPLUS<U,N>(val.getTTypeNameHV());	
	return TTypeName<U,N>(pHV);
}

// POWER

template <typename U, int N>
struct POW : public UnTTypeNameHV<U,N>
{
	POW(const U& val, TTypeNameHV<U,N>* pOp):UnTTypeNameHV<U,N>(val,pOp){}
	POW(TTypeNameHV<U,N>* pOp):UnTTypeNameHV<U,N>(pOp){}
	unsigned int eval(const unsigned int k)
	{
		unsigned int l=this->opEval(k);
		for(unsigned int i=this->length();i<l;++i) this->val(i)=this->opVal(i);
		return this->length()=l;
	}
private:
	void operator=(const POW<U,N>&){} // not allowed
};
template <typename U, int N, typename V>
struct POW1 : public UnTTypeNameHV<U,N>
{
	POW1(const U& val, TTypeNameHV<U,N>* pOp2):UnTTypeNameHV<U,N>(val,pOp2){}
	POW1(TTypeNameHV<U,N>* pOp2):UnTTypeNameHV<U,N>(pOp2){}
	unsigned int eval(const unsigned int k)
	{
		unsigned int l=this->opEval(k);
		for(unsigned int i=this->length();i<l;++i) this->val(i)=this->opVal(i);
		return this->length()=l;
	}
private:
	void operator=(const POW1<U,N,V>&){} // not allowed
};
template <typename U, int N, typename V>
struct POW2 : public UnTTypeNameHV<U,N>
{
	POW2(const U& val, TTypeNameHV<U,N>* pOp1):UnTTypeNameHV<U,N>(val,pOp1){}
	POW2(TTypeNameHV<U,N>* pOp1):UnTTypeNameHV<U,N>(pOp1){}
	unsigned int eval(const unsigned int k)
	{
		unsigned int l=this->opEval(k);
		for(unsigned int i=this->length();i<l;++i) this->val(i)=this->opVal(i);
		return this->length()=l;
	}
private:
	void operator=(const POW2<U,N,V>&){} // not allowed
};
template <typename U, int N>
TTypeName<U,N> pow(const TTypeName<U,N>& val1, const TTypeName<U,N>& val2)
{
	TTypeName<U,N> tmp(exp(val2*log(val1)));
	TTypeNameHV<U,N>* pHV=val1.length()>0 && val2.length()>0 ?
		new POW<U,N>(Op<U>::myPow(val1.val(),val2.val()),tmp.getTTypeNameHV()) :
		new POW<U,N>(tmp.getTTypeNameHV());
	return TTypeName<U,N>(pHV);
}
template <typename U, int N, typename V>
TTypeName<U,N> pow(const V& a, const TTypeName<U,N>& val2)
{
	TTypeName<U,N> tmp(exp(val2*Op<V>::myLog(a)));
	TTypeNameHV<U,N>* pHV=val2.length()>0 ?
		new POW1<U,N,V>(Op<U>::myPow(a,val2.val()), tmp.getTTypeNameHV()) :
		new POW1<U,N,V>(tmp.getTTypeNameHV());
	return TTypeName<U,N>(pHV);
}
template <typename U, int N, typename V>
TTypeName<U,N> pow(const TTypeName<U,N>& val1, const V& b)
{
	TTypeName<U,N> tmp(exp(b*log(val1)));
	TTypeNameHV<U,N>* pHV=val1.length()>0 ?
		new POW2<U,N,V>(Op<U>::myPow(val1.val(),b), tmp.getTTypeNameHV()) :
		new POW2<U,N,V>(tmp.getTTypeNameHV());
	return TTypeName<U,N>(pHV);
}

// SQR

template <typename U, int N>
struct SQR : public UnTTypeNameHV<U,N>
{
	SQR(const U& val, TTypeNameHV<U,N>* pOp):UnTTypeNameHV<U,N>(val,pOp){}
	SQR(TTypeNameHV<U,N>* pOp):UnTTypeNameHV<U,N>(pOp){}
	unsigned int eval(const unsigned int k)
	{
		unsigned int l=this->opEval(k);
		if (0==this->length()) { this->val(0)=Op<U>::mySqr(this->opVal(0)); this->length()=1; }
		for(unsigned int i=this->length();i<l;++i)
		{
			this->val(i)=Op<U>::myZero();
			unsigned int m=(i+1)/2;
			for(unsigned int j=0;j<m;++j) Op<U>::myCadd(this->val(i), this->opVal(i-j)*this->opVal(j));
			Op<U>::myCmul(this->val(i), Op<U>::myTwo());
			if (0==i%2) Op<U>::myCadd(this->val(i), Op<U>::mySqr(this->opVal(m)));
		}
		return this->length()=l;
	}
private:
	void operator=(const SQR<U,N>&){} // not allowed
};
template <typename U, int N>
TTypeName<U,N> sqr(const TTypeName<U,N>& val)
{ 
	TTypeNameHV<U,N>* pHV=val.length()>0 ?
		new SQR<U,N>(Op<U>::mySqr(val.val()), val.getTTypeNameHV()) :
		new SQR<U,N>(val.getTTypeNameHV());
	return TTypeName<U,N>(pHV);
}

// SQRT

template <typename U, int N>
struct SQRT : public UnTTypeNameHV<U,N>
{
	SQRT(const U& val, TTypeNameHV<U,N>* pOp):UnTTypeNameHV<U,N>(val,pOp){}
	SQRT(TTypeNameHV<U,N>* pOp):UnTTypeNameHV<U,N>(pOp){}
	unsigned int eval(const unsigned int k)
	{
		unsigned int l=this->opEval(k);
		if (0==this->length()) { this->val(0)=Op<U>::mySqrt(this->opVal(0)); this->length()=1; }
		for(unsigned int i=this->length();i<l;++i)
		{
			this->val(i)=Op<U>::myZero();
			unsigned int m=(i+1)/2;
			for(unsigned int j=1;j<m;++j) Op<U>::myCadd(this->val(i), this->val(i-j)*this->val(j));
			Op<U>::myCmul(this->val(i), Op<U>::myTwo());
			if (0==i%2) Op<U>::myCadd(this->val(i), Op<U>::mySqr(this->val(m)));
			this->val(i)=(this->opVal(i)-this->val(i))/(Op<U>::myTwo()*this->val(0));
		}
		return this->length()=l;
	}
private:
	void operator=(const SQRT<U,N>&){} // not allowed
};
template <typename U, int N>
TTypeName<U,N> sqrt(const TTypeName<U,N>& val)
{
	TTypeNameHV<U,N>* pHV=val.length()>0 ?
		new SQRT<U,N>(Op<U>::mySqrt(val.val()), val.getTTypeNameHV()):
		new SQRT<U,N>(val.getTTypeNameHV());
	return TTypeName<U,N>(pHV);
}

// EXP

template <typename U, int N>
struct EXP : public UnTTypeNameHV<U,N>
{
	EXP(const U& val, TTypeNameHV<U,N>* pOp):UnTTypeNameHV<U,N>(val,pOp){}
	EXP(TTypeNameHV<U,N>* pOp):UnTTypeNameHV<U,N>(pOp){}
	unsigned int eval(const unsigned int k)
	{
		unsigned int l=this->opEval(k);
		if (0==this->length()) { this->val(0)=Op<U>::myExp(this->opVal(0)); this->length()=1; }
		for(unsigned int i=this->length();i<l;++i)
		{
			this->val(i)=Op<U>::myZero();
			for(unsigned int j=0;j<i;++j)
				Op<U>::myCadd(this->val(i), (Op<U>::myOne()-Op<U>::myInteger(j) / 
					Op<U>::myInteger(i))*this->opVal(i-j)*this->val(j));
		}
		return this->length()=l;
	}
private:
	void operator=(const EXP<U,N>&){} // not allowed
};
template <typename U, int N>
TTypeName<U,N> exp(const TTypeName<U,N>& val)
{
	TTypeNameHV<U,N>* pHV=val.length()>0 ?
		new EXP<U,N>(Op<U>::myExp(val.val()), val.getTTypeNameHV()):
		new EXP<U,N>(val.getTTypeNameHV());
	return TTypeName<U,N>(pHV);
}

// LOG

template <typename U, int N>
struct LOG : public UnTTypeNameHV<U,N>
{
	LOG(const U& val, TTypeNameHV<U,N>* pOp):UnTTypeNameHV<U,N>(val,pOp){}
	LOG(TTypeNameHV<U,N>* pOp):UnTTypeNameHV<U,N>(pOp){}
	unsigned int eval(const unsigned int k)
	{
		unsigned int l=this->opEval(k);
		if (0==this->length()) { this->val(0)=Op<U>::myLog(this->opVal(0)); this->length()=1; }
		for(unsigned int i=this->length();i<l;++i)
		{
			this->val(i)=this->opVal(i);
			for(unsigned int j=1;j<i;++j)
				Op<U>::myCsub(this->val(i), (Op<U>::myOne()-Op<U>::myInteger(j) /
					Op<U>::myInteger(i))*this->opVal(j)*this->val(i-j));
			Op<U>::myCdiv(this->val(i), this->opVal(0));
		}
		return this->length()=l;
	}
private:
	void operator=(const LOG<U,N>&){} // not allowed
};
template <typename U, int N>
TTypeName<U,N> log(const TTypeName<U,N>& val)
{
	TTypeNameHV<U,N>* pHV=val.length()>0 ?
		new LOG<U,N>(Op<U>::myLog(val.val()), val.getTTypeNameHV()):
		new LOG<U,N>(val.getTTypeNameHV());
	return TTypeName<U,N>(pHV);
}

// SIN

template <typename U, int N>
struct SIN : public UnTTypeNameHV<U,N>
{
	U m_COS[N];
	SIN(const U& val, TTypeNameHV<U,N>* pOp):UnTTypeNameHV<U,N>(val,pOp){m_COS[0]=Op<U>::myCos(this->opVal(0));}
	SIN(TTypeNameHV<U,N>* pOp):UnTTypeNameHV<U,N>(pOp){}
	unsigned int eval(const unsigned int k)
	{
		unsigned int l=this->opEval(k);
		if (0==this->length()) 
		{ 
			this->val(0)=Op<U>::mySin(this->opVal(0));
			m_COS[0]=Op<U>::myCos(this->opVal(0));
			this->length()=1; 
		}
		for(unsigned int i=this->length();i<l;++i)
		{
			this->val(i)=Op<U>::myZero();
			for(unsigned int j=0;j<i;++j)
				Op<U>::myCadd(this->val(i), Op<U>::myInteger(j+1)*m_COS[i-1-j]*this->opVal(j+1));
			Op<U>::myCdiv(this->val(i), Op<U>::myInteger(i));
			m_COS[i]=Op<U>::myZero();
			for(unsigned int j=0;j<i;++j)
				Op<U>::myCsub(m_COS[i], Op<U>::myInteger(j+1)*this->val(i-1-j)*this->opVal(j+1));
			Op<U>::myCdiv(m_COS[i], Op<U>::myInteger(i));
		}
		return this->length()=l;
	}
private:
	void operator=(const SIN<U,N>&){} // not allowed
};
template <typename U, int N>
TTypeName<U,N> sin(const TTypeName<U,N>& val)
{
	TTypeNameHV<U,N>* pHV=val.length()>0 ?
		new SIN<U,N>(Op<U>::mySin(val.val()), val.getTTypeNameHV()):
		new SIN<U,N>(val.getTTypeNameHV());
	return TTypeName<U,N>(pHV);
}

// COS

template <typename U, int N>
struct COS : public UnTTypeNameHV<U,N>
{
	U m_SIN[N];
	COS(const U& val, TTypeNameHV<U,N>* pOp):UnTTypeNameHV<U,N>(val,pOp){m_SIN[0]=Op<U>::mySin(this->opVal(0));}
	COS(TTypeNameHV<U,N>* pOp):UnTTypeNameHV<U,N>(pOp){}
	unsigned int eval(const unsigned int k)
	{
		unsigned int l=this->opEval(k);
		if (0==this->length()) 
		{ 
			this->val(0)=Op<U>::myCos(this->opVal(0));
			m_SIN[0]=Op<U>::mySin(this->opVal(0));
			this->length()=1; 
		}
		for(unsigned int i=this->length();i<l;++i)
		{
			this->val(i)=Op<U>::myZero();
			for(unsigned int j=0;j<i;++j)
				Op<U>::myCsub(this->val(i), Op<U>::myInteger(j+1)*m_SIN[i-1-j]*this->opVal(j+1));
			Op<U>::myCdiv(this->val(i), Op<U>::myInteger(i));
			m_SIN[i]=Op<U>::myZero();
			for(unsigned int j=0;j<i;++j)
				Op<U>::myCadd(m_SIN[i], Op<U>::myInteger(j+1)*this->val(i-1-j)*this->opVal(j+1));
			Op<U>::myCdiv(m_SIN[i], Op<U>::myInteger(i));
		}
		return this->length()=l;
	}
private:
	void operator=(const COS<U,N>&){} // not allowed
};
template <typename U, int N>
TTypeName<U,N> cos(const TTypeName<U,N>& val)
{
	TTypeNameHV<U,N>* pHV=val.length()>0 ?
		new COS<U,N>(Op<U>::myCos(val.val()), val.getTTypeNameHV()):
		new COS<U,N>(val.getTTypeNameHV());
	return TTypeName<U,N>(pHV);
}

// TAN

template <typename U, int N>
struct TAN : public BinTTypeNameHV<U,N>
{
	TAN(const U& val, TTypeNameHV<U,N>* pOp, TTypeNameHV<U,N>* pSqrCos):BinTTypeNameHV<U,N>(val,pOp,pSqrCos){}
	TAN(TTypeNameHV<U,N>* pOp, TTypeNameHV<U,N>* pSqrCos):BinTTypeNameHV<U,N>(pOp,pSqrCos){}
	unsigned int eval(const unsigned int k)
	{
		unsigned int l=std::min(this->op1Eval(k),this->op2Eval(k));
		if (0==this->length()) { this->val(0)=Op<U>::myTan(this->op1Val(0)); this->length()=1; }
		for(unsigned int i=this->length();i<l;++i)
		{
			this->val(i)=Op<U>::myZero();
			for(unsigned int j=1;j<i;++j)
				Op<U>::myCadd(this->val(i), Op<U>::myInteger(j)*this->val(j)*this->op2Val(i-j));
			this->val(i)=(this->op1Val(i)-this->val(i)/Op<U>::myInteger(i))/this->op2Val(0);
		}
		return this->length()=l;
	}
private:
	void operator=(const TAN<U,N>&){} // not allowed
};
template <typename U, int N>
TTypeName<U,N> tan(const TTypeName<U,N>& val)
{ 
	TTypeName<U,N> tmp(sqr(cos(val)));
	TTypeNameHV<U,N>* pHV=val.length()>0 ?
		new TAN<U,N>(tan(val.val()), val.getTTypeNameHV(), tmp.getTTypeNameHV()):
		new TAN<U,N>(val.getTTypeNameHV(), tmp.getTTypeNameHV());
	return TTypeName<U,N>(pHV);
}

// ASIN

template <typename U, int N>
struct ASIN : public BinTTypeNameHV<U,N>
{
	ASIN(const U& val, TTypeNameHV<U,N>* pOp, TTypeNameHV<U,N>* pSqrt):BinTTypeNameHV<U,N>(val,pOp,pSqrt){}
	ASIN(TTypeNameHV<U,N>* pOp, TTypeNameHV<U,N>* pSqrt):BinTTypeNameHV<U,N>(pOp,pSqrt){}
	unsigned int eval(const unsigned int k)
	{
		unsigned int l=std::min(this->op1Eval(k),this->op2Eval(k));
		if (0==this->length()) { this->val(0)=Op<U>::myAsin(this->op1Val(0)); this->length()=1; }
		for(unsigned int i=this->length();i<l;++i)
		{
			this->val(i)=Op<U>::myZero();
			for(unsigned int j=1;j<i;++j)
				Op<U>::myCadd(this->val(i), Op<U>::myInteger(j)*this->val(j)*this->op2Val(i-j));
			this->val(i)=(this->op1Val(i)-this->val(i)/Op<U>::myInteger(i))/this->op2Val(0);
		}
		return this->length()=l;
	}
private:
	void operator=(const ASIN<U,N>&){} // not allowed
};
template <typename U, int N>
TTypeName<U,N> asin(const TTypeName<U,N>& val)
{
	TTypeName<U,N> tmp(sqrt(Op<U>::myOne()-sqr(val)));
	TTypeNameHV<U,N>* pHV=val.length()>0 ?
		new ASIN<U,N>(Op<U>::myAsin(val.val()), val.getTTypeNameHV(), tmp.getTTypeNameHV()):
		new ASIN<U,N>(val.getTTypeNameHV(), tmp.getTTypeNameHV());
	return TTypeName<U,N>(pHV);
}

// ACOS

template <typename U, int N>
struct ACOS : public BinTTypeNameHV<U,N>
{
	ACOS(const U& val, TTypeNameHV<U,N>* pOp, TTypeNameHV<U,N>* pSqrt):BinTTypeNameHV<U,N>(val,pOp,pSqrt){}
	ACOS(TTypeNameHV<U,N>* pOp, TTypeNameHV<U,N>* pSqrt):BinTTypeNameHV<U,N>(pOp,pSqrt){}
	unsigned int eval(const unsigned int k)
	{
		unsigned int l=std::min(this->op1Eval(k),this->op2Eval(k));
		if (0==this->length()) { this->val(0)=Op<U>::myAcos(this->op1Val(0)); this->length()=1; }
		for(unsigned int i=this->length();i<l;++i)
		{
			this->val(i)=Op<U>::myZero();
			for(unsigned int j=1;j<i;++j)
				Op<U>::myCadd(this->val(i), Op<U>::myInteger(j)*this->val(j)*this->op2Val(i-j));
			this->val(i)=Op<U>::myNeg((this->op1Val(i)+this->val(i)/Op<U>::myInteger(i))/this->op2Val(0));
		}
		return this->length()=l;
	}
private:
	void operator=(const ACOS<U,N>&){} // not allowed
};
template <typename U, int N>
TTypeName<U,N> acos(const TTypeName<U,N>& val)
{
	TTypeName<U,N> tmp(sqrt(Op<U>::myOne()-sqr(val)));
	TTypeNameHV<U,N>* pHV=val.length()>0 ?
		new ACOS<U,N>(Op<U>::myAcos(val.val()), val.getTTypeNameHV(), tmp.getTTypeNameHV()):
		new ACOS<U,N>(val.getTTypeNameHV(), tmp.getTTypeNameHV());
	return TTypeName<U,N>(pHV);
}

// ATAN

template <typename U, int N>
struct ATAN : public BinTTypeNameHV<U,N>
{
	ATAN(const U& val, TTypeNameHV<U,N>* pOp, TTypeNameHV<U,N>* p1pSqr):BinTTypeNameHV<U,N>(val,pOp,p1pSqr){}
	ATAN(TTypeNameHV<U,N>* pOp, TTypeNameHV<U,N>* p1pSqr):BinTTypeNameHV<U,N>(pOp,p1pSqr){}
	unsigned int eval(const unsigned int k)
	{
		unsigned int l=std::min(this->op1Eval(k),this->op2Eval(k));
		if (0==this->length()) { this->val(0)=Op<U>::myAtan(this->op1Val(0)); this->length()=1; }
		for(unsigned int i=this->length();i<l;++i)
		{
			this->val(i)=Op<U>::myZero();
			for(unsigned int j=1;j<i;++j)
				Op<U>::myCadd(this->val(i), Op<U>::myInteger(j)*this->val(j)*this->op2Val(i-j));
			this->val(i)=(this->op1Val(i)-this->val(i)/Op<U>::myInteger(i))/this->op2Val(0);
		}
		return this->length()=l;
	}
private:
	void operator=(const ATAN<U,N>&){} // not allowed
};
template <typename U, int N>
TTypeName<U,N> atan(const TTypeName<U,N>& val)
{ 
	TTypeName<U,N> tmp(Op<U>::myOne()+sqr(val));
	TTypeNameHV<U,N>* pHV=val.length()>0 ?
		new ATAN<U,N>(Op<U>::myAtan(val.val()), val.getTTypeNameHV(), tmp.getTTypeNameHV()):
		new ATAN<U,N>(val.getTTypeNameHV(), tmp.getTTypeNameHV());
	return TTypeName<U,N>(pHV);
}

// Ned's diff operator

template <typename U, int N>
struct DIFF : public UnTTypeNameHV<U,N>
{
	int m_b;
	DIFF(const U& val, TTypeNameHV<U,N>* pOp, const int b):UnTTypeNameHV<U,N>(val,pOp),m_b(b){}
	DIFF(TTypeNameHV<U,N>* pOp, const int b):UnTTypeNameHV<U,N>(pOp),m_b(b){}
	unsigned int eval(const unsigned int k)
	{
// IN ORDER TO COMPUTE i'th ORDER COEFFICIENTS OF diff(m_o1,b)
// WE NEED (i+b)'th ORDER COEFFICIENTS OF TaylorOp.
		unsigned int l=this->opEval(k+m_b);
// NOW WE SHOULD HAVE op1 EXPANDED TO DEGREE (i+b)'th ORDER,
// WE CAN PROCEED TO COMPUTE UP TO i'TH ORDER COEFFICIENTS
// OF diff(op1,b).
		if (this->length()+m_b<l)
		{
			for(unsigned int i=this->length();i<l-m_b;++i)
			{
				unsigned int fact=1;
				for(unsigned int j=i+m_b;j>i;--j){ fact*=j; }
				this->val(i)=this->opVal(i+m_b)*fact;
			}
			this->length()=l-m_b;
		}
		return this->length();
	}
private:
	void operator=(const DIFF<U,N>&){} // not allowed
};
template <typename U, int N>
TTypeName<U,N> diff(const TTypeName<U,N>& val, const int b)
{
// IF THE ARGUMENT TAYLOR EXPANSION TaylorOp HAS BEEN EVALUATED TO
// DEGREE i THEN WE CAN EVALUATE THE ZERO ORDER VALUE OF
// diff(TaylorOp,i).
// THIS FUNCTION EVALUATES THE 0.ORDER COEFFICIENT
	TTypeNameHV<U,N>* pHV=0;
	if (val.length()>b)
	{
		unsigned int fact=1;
		for(unsigned int j=b;j>1;--j){ fact*=j; }
		pHV=new DIFF<U,N>(val[b]*fact, val.getTTypeNameHV(), b);
	}
	else
	{
		pHV=new DIFF<U,N>(val.getTTypeNameHV(), b);
	}
	return TTypeName<U,N>(pHV);
}


template <typename U, int N> struct Op< TTypeName<U,N> >
{
	typedef TTypeName<U,N> V;
	typedef typename Op<U>::Base Base;
	static Base myInteger(const int i) { return Base(i); }
	static Base myZero() { return myInteger(0); }
	static Base myOne() { return myInteger(1);}
	static Base myTwo() { return myInteger(2); }
	static V myPos(const V& x) { return +x; }
	static V myNeg(const V& x) { return -x; }
	static V& myCadd(V& x, const V& y) { return x+=y; }
	static V& myCsub(V& x, const V& y) { return x-=y; }
	static V& myCmul(V& x, const V& y) { return x*=y; }
	static V& myCdiv(V& x, const V& y) { return x/=y; }
	static V myInv(const V& x) { return myOne()/x; }
	static V mySqr(const V& x) { return x*x; }
	static V myPow(const V& x, const int n) { return ::pow(x,n); }
	static V myPow(const V& x, const V& y) { return ::pow(x,y); }
	static V mySqrt(const V& x) { return ::sqrt(x); }
	static V myLog(const V& x) { return ::log(x); }
	static V myExp(const V& x) { return ::exp(x); }
	static V mySin(const V& x) { return ::sin(x); }
	static V myCos(const V& x) { return ::cos(x); }
	static V myTan(const V& x) { return ::tan(x); }
	static V myAsin(const V& x) { return ::asin(x); }
	static V myAcos(const V& x) { return ::acos(x); }
	static V myAtan(const V& x) { return ::atan(x); }
	static bool myEq(const V& x, const V& y) { return x==y; }
	static bool myNe(const V& x, const V& y) { return x!=y; }
	static bool myLt(const V& x, const V& y) { return x<y; }
	static bool myLe(const V& x, const V& y) { return x<=y; }
	static bool myGt(const V& x, const V& y) { return x>y; }
	static bool myGe(const V& x, const V& y) { return x>=y; }
};


#endif
