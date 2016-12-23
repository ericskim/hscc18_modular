/*363:*/
#line 390 "./fadbad.w"

#ifndef FADBADAD_H
#define FADBADAD_H
#include "fadbad_ad.h"
#include "fadbad_advar.h"
#include "allad.h"

namespace vnodelp{
/*362:*/
#line 372 "./fadbad.w"

class FADBAD_AD:public AD
{
public:
FADBAD_AD(int n,Tfunction f,TFfunction tf)
:AD(n,new FadbadODE(n,f),new FadbadVarODE(n,tf)),max_order(MaxLength-2){}

FADBAD_AD(int n,Tfunction f,TFfunction tf,void*p)
:AD(n,new FadbadODE(n,f,p),new FadbadVarODE(n,tf,p)),
max_order(MaxLength-2){}

virtual int getMaxOrder()const{return max_order;}

private:
const int max_order;
};

/*:362*/
#line 398 "./fadbad.w"

}
#endif
#line 123 "./vnode.w"

/*:363*/
