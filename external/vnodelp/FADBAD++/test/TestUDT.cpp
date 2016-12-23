
// USER DEFINED TYPE TEST:

#include "fadbad.h"
#include <iostream>

#include "TestUDT.h"

class UDT
{
public:
	UDT(){/*std::cout<<"UDT()"<<std::endl;*/} // default constructor
	UDT(const UDT&){/*std::cout<<"UDT(UDT&)"<<std::endl;*/} // copy constructor
	UDT(const double&){/*std::cout<<"UDT(UDT&)"<<std::endl;*/} // copy constructor
	UDT(const int){/*std::cout<<"UDT(int)"<<std::endl;*/} // initialize-with-int constructor
	UDT& operator=(const UDT&){ return *this; }
	UDT& operator=(const double&){ return *this; }
};

UDT operator + (const UDT&, const UDT&) {/* std::cout<<"UDT+UDT"<<std::endl;*/return UDT(); }
UDT operator - (const UDT&, const UDT&) {/* std::cout<<"UDT-UDT"<<std::endl;*/return UDT(); }
UDT operator * (const UDT&, const UDT&) {/* std::cout<<"UDT*UDT"<<std::endl;*/return UDT(); }
UDT operator / (const UDT&, const UDT&) {/* std::cout<<"UDT/UDT"<<std::endl;*/return UDT(); }
UDT Exp(const UDT&) {/* std::cout<<"Exp(UDT)"<<std::endl;*/return UDT(); }
UDT Sqr(const UDT&) {/* std::cout<<"Sqr(UDT)"<<std::endl;*/return UDT(); }
UDT Sin(const UDT&) {/* std::cout<<"Sin(UDT)"<<std::endl;*/return UDT(); }
UDT Cos(const UDT&) {/* std::cout<<"Cos(UDT)"<<std::endl;*/return UDT(); }

UDT operator - (const UDT&) {/* std::cout<<"-UDT"<<std::endl;*/return UDT(); }


template <>
struct Op<UDT>
{
	typedef UDT Base;
	static Base myInteger(const int i) { return Base(i); }
	static Base myZero() { return Base(0); }
	static Base myOne() { return Base(1);}
	static Base myTwo() { return Base(2); }
	static UDT myPos(const UDT& /*x*/) { return UDT(); }
	static UDT myNeg(const UDT& /*x*/) { return UDT(); }
	static UDT& myCadd(UDT& x, const UDT& y) { x=x+y; return x; }
	static UDT& myCsub(UDT& x, const UDT& y) { x=x-y; return x; }
	static UDT& myCmul(UDT& x, const UDT& y) { x=x*y; return x; }
	static UDT& myCdiv(UDT& x, const UDT& y) { x=x/y; return x; }
	static UDT myInv(const UDT& x) { return myOne()/x; }
	static UDT mySqr(const UDT& x) { return ::Sqr(x); }
	static UDT myPow(const UDT& /*x*/, const int /*n*/) { throw; return UDT(); }
	static UDT myPow(const UDT& /*x*/, const UDT& /*y*/) { throw; return UDT(); }
	static UDT mySqrt(const UDT& /*x*/) { throw; return UDT(); }
	static UDT myLog(const UDT& /*x*/) { throw; return UDT(); }
	static UDT myExp(const UDT& x) { return ::Exp(x); }
	static UDT mySin(const UDT& x) { return ::Sin(x); }
	static UDT myCos(const UDT& x) { return ::Cos(x); }
	static UDT myTan(const UDT& /*x*/) { throw; return UDT(); }
	static UDT myAsin(const UDT& /*x*/) { throw; return UDT(); }
	static UDT myAcos(const UDT& /*x*/) { throw; return UDT(); }
	static UDT myAtan(const UDT& /*x*/) { throw; return UDT(); }
	static bool myEq(const UDT& /*x*/, const UDT& /*y*/) { throw; return false; }
	static bool myNe(const UDT& /*x*/, const UDT& /*y*/) { throw; return false; }
	static bool myLt(const UDT& /*x*/, const UDT& /*y*/) { throw; return false; }
	static bool myLe(const UDT& /*x*/, const UDT& /*y*/) { throw; return false; }
	static bool myGt(const UDT& /*x*/, const UDT& /*y*/) { throw; return false; }
	static bool myGe(const UDT& /*x*/, const UDT& /*y*/) { throw; return false; }
};

#include "badiff.h" // Backwards automatic differentiation
#include "fadiff.h" // Forwards automatic differentiation
#include "ffadiff.h" // Fast-Forward automatic differentiation
#include "tadiff.h" // Taylor expansion

template< int N, class U>
class RecTestTypes
{
	RecTestTypes< N-1,U > a;
	RecTestTypes< N-1,F<U> > b;
	RecTestTypes< N-1,FF<U,1> > c;
	RecTestTypes< N-1,B<U> > d;
	RecTestTypes< N-1,T<U> > e;
};

template< class U >
class RecTestTypes<0,U>
{
public:
	RecTestTypes()
	{
		U a;
		double x=1.23;
		typename Op<U>::Base t=Op<U>::myTwo();
		// Ensure we can use some minimal set of operations:
		F<U> b(a);
		F<U> b1(b);
		F<U> b2=x;
		FF<U,1> c(a);
		FF<U,1> c1(c);
		FF<U,1> c2=x;
//		B<U> d(a);
//		B<U> d1(d);
//		B<U> d2(x);
		T<U> d(a);
		T<U> d1(d);
		T<U> d2=x;
		T<U> e(a);
		T<U> e1(e);
		T<U> e2=x;
		a+a;b+b;c+c;d+d;e+e;
//		b+a;c+a;d+a;e+a;
//		a+b;a+c;a+d;a+e;
		b+t;c+t;d+t;e+t;
		t+b;t+c;t+d;t+e;
		b+x;c+x;d+x;e+x;
		x+b;x+c;x+d;x+e;
		a-a;b-b;c-c;d-d;e-e;
//		b-a;c-a;d-a;e-a;
//		a-b;a-c;a-d;a-e;
		b-t;c-t;d-t;e-t;
		t-b;t-c;t-d;t-e;
		b-x;c-x;d-x;e-x;
		x-b;x-c;x-d;x-e;
		a*a;b*b;c*c;d*d;e*e;
//		b*a;c*a;d*a;e*a;
//		a*b;a*c;a*d;a*e;
		b*t;c*t;d*t;e*t;
		t*b;t*c;t*d;t*e;
		b*x;c*x;d*x;e*x;
		x*b;x*c;x*d;x*e;
		a/a;b/b;c/c;d/d;e/e;
//		b/a;c/a;d/a;e/a;
//		a/b;a/c;a/d;a/e;
		b/t;c/t;d/t;e/t;
		t/b;t/c;t/d;t/e;
		b/x;c/x;d/x;e/x;
		x/b;x/c;x/d;x/e;
		a=a;b=b;c=c;d=d;e=e;
		b=a;c=a;d=a;e=a;
		b=t;c=t;d=t;e=t;
		b=x;c=x;d=x;e=x;
	}
};

void TestUDT::run(IReportLog& rlog)
{
	RecTestTypes<0,UDT> t0;
	rlog.succeeded() << "User defined type test level 0 succeeded" << std::endl;
	RecTestTypes<1,UDT> t1;
	rlog.succeeded() << "User defined type test level 1 succeeded" << std::endl;
//	RecTestTypes<2,UDT> t2;
//	rlog.succeeded() << "User defined type test level 2 succeeded" << std::endl;
}	



