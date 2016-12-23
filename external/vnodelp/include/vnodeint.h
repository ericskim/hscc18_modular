/*333:*/
#line 913 "./integ.w"

#ifndef VNODEINT_H
#define VNODEINT_H

#include <algorithm> 
#include "miscfuns.h"
#include "vnodeinterval.h"
#include "vnoderound.h"
#include "vector_matrix.h"
#include "ad_ode.h"
#include "ad_var.h"
#include "allad.h"
#include "solution.h"
#include "control.h"
#include "matrixinverse.h"
#include "iho.h"
#include "hoe.h"
#include "vtiming.h"
#include "debug.h"



namespace vnodelp{
/*295:*/
#line 6 "./integ.w"


typedef enum{on,off}stepAction;



class VNODE
{
public:
int steps;
VNODE(AD*ad);

void integrate(interval&t0,iVector&y0,const interval&t_end);

/*331:*/
#line 863 "./integ.w"



void setTols(double a,double r= 0)
{


control->atol= a;
control->rtol= r;
}



void setOrder(int p)
{
control->order= p;
}


void setHmin(double h)
{
control->hmin= h;
}




void setOneStep(stepAction action)
{
if(action==on)
control->interrupt= before_accept;
else
control->interrupt= no;

}

void setFirstEntry(){
control->ind= first_entry;
}








/*:331*/
#line 20 "./integ.w"

/*329:*/
#line 808 "./integ.w"


unsigned int getMaxOrder()const{return ad->getMaxOrder();}

bool successful()const
{
if(control->ind==success||control->ind==first_entry)
return true;

return false;
}

/*:329*//*330:*/
#line 821 "./integ.w"


double getStepsize()const
{
return midpoint(h_accepted);
}

const iVector&getAprioriEncl()const
{
return hoe->getApriori();
}



const interval&getT()const
{
return hoe->getT();
}


double getGlobalExcess(){

width(tp,iho->getGlobalExcess());
return inf_normV(tp);
}



double getGlobalExcess(int i){


width(tp,iho->getGlobalExcess());
return tp[i];
}

int getNoSteps()const{

return steps;
}

/*:330*/
#line 21 "./integ.w"


~VNODE();

private:
void acceptSolution(interval&t0,iVector&y0);
double compHstart(const interval&t0,const iVector&y0);

private:
int direction;
interval t_trial,Tj,h_accepted;
iVector temp;
double h_start,h_min;
pVector tp;

Control*control;
HOE*hoe;
IHO*iho;

AD*ad;
};



/*:295*/
#line 936 "./integ.w"

}

#endif
/*:333*/
