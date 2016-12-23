/*227:*/
#line 842 "./hoe.w"

#include <cmath> 
#include <cassert> 
#include <algorithm> 

#include "vnodeinterval.h"
#include "vnoderound.h"
#include "basiclinalg.h"
#include "intvfuncs.h"
#include "control.h"
#include "solution.h"
#include "allad.h"
#include "hoe.h"

using namespace v_blas;



namespace vnodelp{

/*219:*/
#line 678 "./hoe.w"

HOE::HOE(int n):one(interval(0,1))
{
sizeV(term,n);
sizeV(p,n);
sizeV(u,n);
sizeV(v,n);

apriori_trial= new Apriori(n);
apriori= new Apriori(n);
assert(apriori&&apriori_trial);
tayl_coeff= 0;
control= 0;
}


HOE::~HOE()
{
delete apriori;
delete apriori_trial;
}



/*:219*/
#line 862 "./hoe.w"

/*217:*/
#line 630 "./hoe.w"

void HOE::compAprioriEnclosure(const interval&t0,const iVector&y0,
bool&info)
{


















if(fabs(h_trial)<=control->hmin){
info= false;
return;
}
/*206:*/
#line 184 "./hoe.w"



















tayl_coeff->set(t0,y0,h_trial,order_trial-1);
tayl_coeff->compTerms();
interval t_enc= (t0-t0)/h_trial+one;
tayl_coeff->getTerm(p,order_trial-1);
for(int i= order_trial-2;i>=0;i--)
{
scaleV(p,t_enc);
tayl_coeff->getTerm(term,i);

addViVi(p,term);
}


/*:206*/
#line 656 "./hoe.w"
;
/*207:*/
#line 235 "./hoe.w"















round_up();
double tol= inf_normV(y0)*control->rtol+control->atol;
setV(u,h_trial*interval(-tol/2,tol/2));
addViVi(apriori_trial->y,p,u);


assert(interior(y0,apriori_trial->y));

/*:207*/
#line 657 "./hoe.w"
;
/*210:*/
#line 301 "./hoe.w"
















tayl_coeff->set(t0+one*h_trial,apriori_trial->y,h_trial,order_trial);
tayl_coeff->compTerms();
tayl_coeff->getTerm(v,order_trial);




/*:210*//*212:*/
#line 369 "./hoe.w"











interval beta= comp_beta(v,u,order_trial);
assert(inf(beta)>=0);
if(inf(beta)<1){

if(h_trial> 0)
round_down();
else
round_up();
h= inf(beta)*h_trial;
}
else
h= h_trial;



/*:212*//*214:*/
#line 444 "./hoe.w"

if(inf(t0)<sup(t0))
{

























interval tt= t0-t0;

while(fabs(h)> control->hmin)
{
t_enc= (tt+one*h)/h_trial;
scaleV(v,pow(t_enc,order_trial));
if(subseteq(v,u))
break;
h= 0.9*h;

}
}
/*:214*/
#line 658 "./hoe.w"
;
if(fabs(h)<=control->hmin)
{info= false;
return;
}
/*215:*/
#line 558 "./hoe.w"










double td;
if(h> 0)
{
td= inf(t0);
round_down();
t_trial= td+h;
apriori_trial->t= interval(td,t_trial);
}
else
{
td= sup(t0);
round_up();
t_trial= td+h;
apriori_trial->t= interval(t_trial,td);
}
assert(subseteq(t0,apriori_trial->t));







/*:215*/
#line 663 "./hoe.w"
;
/*216:*/
#line 607 "./hoe.w"









h_next= inf(beta)*h_trial;




/*:216*/
#line 664 "./hoe.w"
;
info= true;
}







/*:217*/
#line 863 "./hoe.w"

/*220:*/
#line 704 "./hoe.w"

void HOE::acceptSolution()
{











apriori->t= apriori_trial->t;
assignV(apriori->y,apriori_trial->y);
h_trial= h_next;
}




/*:220*/
#line 864 "./hoe.w"

/*224:*/
#line 806 "./hoe.w"








interval HOE::comp_beta(const iVector&v,const iVector&u,int k)
{
double gamma= v_blas::compH(v,u);

interval i_gamma= gamma;

interval i_pw= interval(1.0)/interval(double(k));


interval beta= pow(i_gamma,i_pw);

return beta;
}

/*:224*/
#line 865 "./hoe.w"

}



#line 1 "./iho.w"
/*:227*/
