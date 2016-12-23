/*334:*/
#line 943 "./integ.w"

#include <cmath> 
#include <algorithm> 
#include <ostream> 
using namespace std;

#include "vnodeinterval.h"
#include "vector_matrix.h"
#include "matrixinverse.h"
#include "basiclinalg.h"
#include "vnodeint.h"


namespace vnodelp{

/*327:*/
#line 780 "./integ.w"


VNODE::VNODE(AD*a)
{
steps= 0;
direction= 0;
int n= a->size;
sizeV(temp,n);
sizeV(tp,n);
hoe= new HOE(n);
iho= new IHO(n);
control= new Control;
assert(a);
ad= a;
}
/*:327*/
#line 958 "./integ.w"

/*328:*/
#line 796 "./integ.w"

VNODE::~VNODE()
{
delete control;
delete iho;
delete hoe;
}


/*:328*/
#line 959 "./integ.w"

/*314:*/
#line 353 "./integ.w"

double compHmin(const interval&t0,const interval&t_end)
{

















double t= std::max(v_bias::mag(t0),v_bias::mag(t_end));
double t_next= nextafter(t,t+1);
v_bias::round_up();
t= t_next-t;
t= std::max(t,v_bias::width(t_end));
return t;
}

/*:314*/
#line 960 "./integ.w"

/*315:*/
#line 392 "./integ.w"

double VNODE::compHstart(const interval&t0,const iVector&y0)
{



















int k= control->order;
ad->tayl_coeff_ode->set(t0,y0,1,k+1);
ad->tayl_coeff_ode->compTerms();

ad->tayl_coeff_ode->getTerm(temp,k+1);

double tol= control->rtol*inf_normV(y0)+control->atol;
double t= (k+1)*inf_normV(temp);

/*316:*/
#line 436 "./integ.w"

if(t==0)
t= std::max(control->atol,control->rtol);



/*:316*/
#line 423 "./integ.w"


t= tol/t;
t= std::pow(t,1.0/k);

return t;
}


/*:315*/
#line 961 "./integ.w"

/*296:*/
#line 64 "./integ.w"

void VNODE::integrate(interval&t0,iVector&y0,const interval&t_end)
{

/*298:*/
#line 104 "./integ.w"

if(control->ind==failure)
{
vnodeMessage("Previous call to integrate() failed");
vnodeMessage("Call setFirstEntry() before calling integrate");
return;
}

/*:298*//*299:*/
#line 117 "./integ.w"


if(!v_bias::finite_interval(t0)||!v_bias::finite_interval(t_end))
{
vnodeMessage("t0 and t_end must be finite");
control->ind= failure;
return;
}

for(unsigned int i= 0;i<sizeV(y0);i++)
if(!v_bias::finite_interval(y0[i]))
{
vnodeMessage("y0 must contain finite intervals");
control->ind= failure;
return;
}

/*:299*//*300:*/
#line 136 "./integ.w"

t_trial= t0;
if(t_trial==t_end&&control->ind!=first_entry)
{
vnodeMessage("Set different t_end");
return;
}




/*:300*//*301:*/
#line 148 "./integ.w"

if(control->atol<=0&&control->rtol<=0)
{
vnodeMessage("Set nonnegative tolerances and at least one positive tolerance");
return;
}


/*:301*//*302:*/
#line 157 "./integ.w"

if(control->order<3|control->order> getMaxOrder())
{
vnodeMessage("Set order >= 3 and order <= MAX_ORDER");
control->ind= failure;
return;
}

/*:302*//*303:*/
#line 166 "./integ.w"

if(control->hmin<0)
{
vnodeMessage("Set minimum stepsize >= 0");
return;
}



/*:303*//*304:*/
#line 176 "./integ.w"

if(!v_bias::disjoint(t0,t_end))
{
control->ind= failure;
vnodeMessage("t0 and t_end must be disjoint");
return;
}

/*:304*/
#line 68 "./integ.w"



/*305:*/
#line 191 "./integ.w"

if(control->ind==first_entry)
direction= 0;

int last_direction= direction;


/*:305*//*306:*/
#line 208 "./integ.w"

direction= 1;
if(v_bias::sup(t_end)<v_bias::inf(t0))
direction= -1;


/*:306*//*307:*/
#line 217 "./integ.w"

if(last_direction!=0&&last_direction==-direction)
{
vnodeMessage("Integration direction changed without resetting.\n"
"    Call setFirstEntry() before integrating.");
control->ind= failure;
return;
}


/*:307*/
#line 71 "./integ.w"



if(control->ind==first_entry)
{
/*308:*/
#line 231 "./integ.w"


steps= 0;
/*309:*/
#line 258 "./integ.w"

h_min= control->hmin;
if(control->hmin==0||control->ind==success)
control->hmin= compHmin(t0,t_end);

if(control->hmin> 0&&control->ind==first_entry)
{
double h= compHmin(t0,t_end);
if(control->hmin<h)
{
control->ind= failure;
vnodeMessage("Set larger value for hmin");

return;
}
assert(control->hmin> 0);
}


/*:309*/
#line 234 "./integ.w"

/*310:*/
#line 282 "./integ.w"

if(control->ind==first_entry){
h_start= compHstart(t0,y0);
}

if(h_start<control->hmin)
{
control->ind= failure;
vnodeMessage("Minimum stepsize reached.");
return;
}

if(direction==-1)
{
if(control->ind==first_entry)
h_start= -h_start;
}




/*:310*/
#line 235 "./integ.w"

/*312:*/
#line 325 "./integ.w"


iho->set(control,hoe,ad);

iho->compCoeffs();

iho->init(t0,y0);




/*:312*/
#line 236 "./integ.w"

/*311:*/
#line 308 "./integ.w"


hoe->setTrialStepsize(h_start);

hoe->setTrialOrder(control->order);

hoe->set(control,ad);

hoe->init(t0,y0);


/*:311*/
#line 237 "./integ.w"





/*:308*/
#line 76 "./integ.w"

}

while(t0!=t_end)
{
if(control->ind==first_entry||control->ind==success)
{
/*317:*/
#line 448 "./integ.w"







bool info;
hoe->compAprioriEnclosure(t0,y0,info);
if(info==false)
{
control->ind= failure;
vnodeMessage("Could not validate solution");
return;
}




/*:317*/
#line 83 "./integ.w"

/*318:*/
#line 490 "./integ.w"








Tj= hoe->getTrialT();

if(v_bias::subseteq(t_end,Tj))
{
/*320:*/
#line 533 "./integ.w"

t_trial= t_end;

/*:320*/
#line 502 "./integ.w"

}
else
{
if(v_bias::disjoint(t_end,Tj))
{
/*321:*/
#line 578 "./integ.w"







double theta= v_bias::mag(t_end-t0);

if(v_bias::width(Tj)>=theta/2)
{
double d= theta/2.0;

interval tmp;
if(direction==1)
tmp= interval(0,d);
else
tmp= interval(-d,0);

bool b;
b= v_bias::intersect(Tj,Tj,t0+tmp);
assert(b);

}

if(direction==1)
t_trial= v_bias::sup(Tj);
else
t_trial= v_bias::inf(Tj);



/*:321*/
#line 508 "./integ.w"


}
else
{
interval tmp;
if(direction==1)
tmp= v_bias::sup(Tj);
else
tmp= inf(Tj);

if(v_bias::subseteq(tmp,t_end))
{
/*323:*/
#line 683 "./integ.w"

if(direction==1)
t_trial= v_bias::inf(t_end);
else
t_trial= v_bias::sup(t_end);


/*:323*/
#line 521 "./integ.w"

}
else
assert(0);
}
}

/*:318*/
#line 84 "./integ.w"

}

/*324:*/
#line 693 "./integ.w"



iho->compTightEnclosure(t_trial);


/*:324*/
#line 87 "./integ.w"

/*325:*/
#line 735 "./integ.w"


bool ret= false;
switch(control->ind)
{
case first_entry:
control->ind= success;
case success:

if(control->interrupt==before_accept)
ret= true;
break;
case failure:return;
}

acceptSolution(t0,y0);

if(ret)
return;

/*:325*/
#line 88 "./integ.w"

}


control->hmin= h_min;
}






/*:296*/
#line 962 "./integ.w"

/*326:*/
#line 759 "./integ.w"

void VNODE::acceptSolution(interval&t0,iVector&y0)
{
hoe->acceptSolution();


iho->acceptSolution();

h_accepted= t_trial-t0;

t0= t_trial;
iho->getTightEnclosure(y0);

steps++;
}





/*:326*/
#line 963 "./integ.w"

}






#line 1 "./vnode_int.w"


/*:334*/
