/*196:*/
#line 194 "./adabst.w"

#ifndef ALLAD_H
#define ALLAD_H

#include "ad_ode.h"
#include "ad_var.h"
namespace vnodelp{
/*191:*/
#line 129 "./adabst.w"

class AD
{
public:
AD(int n,AD_ODE*a,AD_VAR*av);
void eval(void*p);

virtual int getMaxOrder()const= 0;

 virtual ~AD(){};
public:
int size;

AD_ODE*tayl_coeff_ode;
AD_VAR*tayl_coeff_var;
}
;

/*:191*/
#line 201 "./adabst.w"

/*192:*/
#line 149 "./adabst.w"

inline AD::AD(int n,AD_ODE
*a,AD_VAR*av):size(n),tayl_coeff_ode(a),tayl_coeff_var(av){}


/*:192*//*193:*/
#line 156 "./adabst.w"

inline void AD::eval(void*p){tayl_coeff_ode->eval(p);tayl_coeff_var->eval(p);}






/*:193*/
#line 202 "./adabst.w"

}
#endif

#line 1 "./vnodecontrol.w"



/*:196*/
