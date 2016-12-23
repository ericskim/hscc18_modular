/*115:*/
#line 6 "./interval_interf.w"


#ifndef VNODEINTERVAL_H
#define VNODEINTERVAL_H

#ifdef PROFIL_VNODE

#include <Interval.h> 
#include <Functions.h> 
#include <LongReal.h> 
#include <LongInterval.h> 

namespace v_bias
{
/*76:*/
#line 13 "./integfun.w"

typedef INTERVAL interval;

/*:76*/
#line 20 "./interval_interf.w"

/*114:*/
#line 5 "./interval_profil.w"


inline double inf(const interval&a){return Inf(a);}

inline double sup(const interval&a){
return Sup(a);}

inline double midpoint(const interval&a){
return Mid(a);
}

inline double width(const interval&b)
{return Diam(b);}

inline double mag(const interval&a){return Abs(a);}


inline bool subseteq(const interval&a,const interval&b)
{return a<=b;}

inline bool interior(const interval&a,const interval&b)
{
return(Inf(a))> (Inf(b))&&((Sup(a))<(Sup(b)));
}

inline bool disjoint(const interval&a,const interval&b)
{
interval c;
return!Intersection(c,a,b);
}

inline bool intersect(interval&c,
const interval&a,const interval&b)
{
return Intersection(c,a,b);
}

inline interval pi(){
return ArcCos(-1.0);
}


inline interval pow(const interval&a,int b){
return Power(a,b);
}

inline interval pow(const interval&a,const interval&b){
return Power(a,b);
}

inline interval exp(const interval&a){
return Exp(a);
}

inline interval log(const interval&a){
return Log(a);
}


inline interval sqr(const interval&a){
return Sqr(a);
}

inline interval sqrt(const interval&a){
return Sqrt(a);
}


inline interval sin(const interval&a){
return Sin(a);
}

inline interval cos(const interval&a){
return Cos(a);
}

inline interval tan(const interval&a){
return Tan(a);
}


inline interval asin(const interval&a){
return ArcSin(a);
}

inline interval acos(const interval&a){
return ArcCos(a);
}

inline interval atan(const interval&a){
return ArcTan(a);
}


inline interval string_to_interval(const char*s)
{
return Enclosure(s);
}








#line 1 "./interval_interf.w"

/*:114*/
#line 21 "./interval_interf.w"

}
#endif


#ifdef FILIB_VNODE

#include <interval/interval.hpp> 
#include <iostream> 
#include <sstream> 
#include <string> 

namespace v_bias{
/*77:*/
#line 17 "./integfun.w"

typedef filib::interval<double> interval;

/*:77*/
#line 34 "./interval_interf.w"

/*113:*/
#line 29 "./interval.w"



inline double inf(const interval&a){return a.inf();}

inline double sup(const interval&a){
return a.sup();}

inline double midpoint(const interval&a){
return a.mid();
}

inline double width(const interval&a)
{return a.diam();}


inline double mag(const interval&a){return a.mag();}


inline bool subseteq(const interval&a,const interval&b)
{return filib::subset(a,b);}

inline bool interior(const interval&a,const interval&b)
{return filib::interior(a,b);}


inline bool disjoint(const interval&a,const interval&b)
{
return filib::disjoint(a,b);
}

inline bool intersect(interval&c,
const interval&a,const interval&b)
{
if(filib::disjoint(a,b))
return false;
c= filib::intersect(a,b);
return true;
}



inline interval pi(){
return filib::interval<double> ::PI();
}



inline interval pow(const interval&a,const interval&b){
return filib::pow(a,b);
}

inline interval pow(const interval&a,int b){
return filib::power(a,b);
}


inline interval exp(const interval&a){
return filib::exp(a);
}

inline interval log(const interval&a){
return filib::log(a);
}

inline interval sqr(const interval&a){
return filib::sqr(a);
}

inline interval sqrt(const interval&a){
return filib::sqrt(a);
}


inline interval sin(const interval&a){
return filib::sin(a);
}

inline interval cos(const interval&a){
return filib::cos(a);
}

inline interval tan(const interval&a){
return filib::tan(a);
}


inline interval asin(const interval&a){
return filib::asin(a);
}

inline interval acos(const interval&a){
return filib::acos(a);
}

inline interval atan(const interval&a){
return filib::atan(a);
}



inline interval string_to_interval(const char*s)
    {
        return interval(s,s);
    }






#line 1 "./interval_profil.w"

/*:113*/
#line 35 "./interval_interf.w"

}
#endif


#endif
#line 1 "./rounding.w"


/*:115*/
