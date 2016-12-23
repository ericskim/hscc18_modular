#ifndef  FADBAD_INTV_INC
#define  FADBAD_INTV_INC

#define MaxLength (MAXORDER+2)



template <> struct Op<v_bias::interval>
{
  
  typedef v_bias::interval Base;

  static Base myInteger(const int i) { return Base(double(i)); }
  static Base myZero() { return myInteger(0); }
  static Base myOne() { return myInteger(1);}
  static Base myTwo() { return myInteger(2); }

  static Base myPos(const Base& x) { return +x; }
  static Base myNeg(const Base& x) { return -x; }
  
  static Base& myCadd(Base& x, const Base& y) { return x+=y; }
  static Base& myCsub(Base& x, const Base& y) { return x-=y; }
  static Base& myCmul(Base& x, const Base& y) { return x*=y; }
  static Base& myCdiv(Base& x, const Base& y) { return x/=y; }
  static Base myInv(const Base& x) { return myOne()/x; }
  static Base mySqr(const Base& x) { return v_bias::sqr(x); }
	
  static Base myPow(const Base& x, const int n) { return v_bias::pow(x,n); }
  static Base myPow(const Base& x, const Base& y) { return v_bias::pow(x,y); }

  static Base mySqrt(const Base& x) { return v_bias::sqrt(x); }
  static Base myLog(const Base& x) { return v_bias::log(x); }
  static Base myExp(const Base& x) { return v_bias::exp(x); }
  static Base mySin(const Base& x) { return v_bias::sin(x); }
  static Base myCos(const Base& x) { return v_bias::cos(x); }
  static Base myTan(const Base& x) { return v_bias::tan(x); }
  static Base myAsin(const Base& x) { return v_bias::asin(x); }
  static Base myAcos(const Base& x) { return v_bias::acos(x); }
  static Base myAtan(const Base& x) { return v_bias::atan(x); }
  
  static bool myEq(const Base& x, const Base& y) { return x==y; }
  static bool myNe(const Base& x, const Base& y) { return x!=y; }
};

#include "tadiff.h"
#include "fadiff.h"



#endif

