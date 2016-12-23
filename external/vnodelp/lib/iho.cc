/*293:*/
#line 1650 "./iho.w"

#include <cmath> 
#include <limits>

#include "vnodeinterval.h"
#include "basiclinalg.h"
#include "solution.h"
#include "control.h"
#include "allad.h"
#include "matrixinverse.h"
#include "iho.h"
#include "debug.h"
#include "miscfuns.h"

namespace vnodelp{
using namespace v_bias;
using namespace v_blas;

extern void sortColumns(pMatrix&B,const pMatrix&A);
extern bool computeQR(pMatrix&B,const pMatrix&A);


/*274:*/
#line 1371 "./iho.w"

IHO::IHO(int n)
{
order_trial= 0;
solution= new Solution(n);

trial_solution= new Solution(n);
matrix_inverse= new MatrixInverse(n);

C_pq= C_qp= 0;
sizeV(y,n);
sizeV(y_pred,n);
sizeV(globalExcess,n);
sizeV(y_pred_point,n);

sizeV(temp,n);
sizeV(temp2,n);
sizeV(x,n);
sizeV(u_next,n);

sizeV(predictor_excess,n);
sizeV(corrector_excess,n);
sizeV(z,n);
sizeV(w,n);
sizeV(gj,n);
sizeV(term,n);
sizeV(d,n);
sizeV(s,n);

sizeM(Fj,n);
sizeM(M,n);
sizeM(Cinv,n);
sizeM(G,n);
sizeM(B,n);
sizeM(S,n);
sizeM(A,n);
sizeM(Q,n);
sizeM(U,n);
sizeM(V,n);
sizeM(Ainv,n);
sizeM(C,n);
sizeM(A_point,n);
}

/*:274*/
#line 1674 "./iho.w"

/*275:*/
#line 1416 "./iho.w"

IHO::~IHO()
{
delete matrix_inverse;
delete trial_solution;
delete solution;
delete[]C_pq;
delete[]C_qp;

}



/*:275*/
#line 1675 "./iho.w"

/*281:*/
#line 1561 "./iho.w"


inline bool isEven(unsigned int k){if(k==0)return true;if(k
&0x01)return false;return true;}



/*:281*/
#line 1676 "./iho.w"

/*279:*/
#line 1504 "./iho.w"


void IHO::compCpq(int p,int q)
{




int tmp= q+p+1;
C_pq[0]= 1.0;

for(int i= 1;i<=p;i++)
{
C_pq[i]= (C_pq[i-1]*double(p-i+1))/double(tmp-i);
}
}




/*:279*/
#line 1677 "./iho.w"

/*280:*/
#line 1541 "./iho.w"

void IHO::compCqp(int p,int q)
{




int tmp= q+p+1;
C_qp[0]= 1.0;
for(int i= 1;i<=q;i++)

C_qp[i]= (-C_qp[i-1]*double(q-i+1))/double(tmp-i);
}



/*:280*/
#line 1678 "./iho.w"

/*282:*/
#line 1578 "./iho.w"






interval IHO::compErrorConstant(int p,int q)
{
interval err_const(1.0);

for(int i= 1;i<=p;i++)
{
err_const*= double(i);
err_const/= (q+i);
}

if(!isEven(q))
err_const= -err_const;

return err_const;
}

/*:282*/
#line 1679 "./iho.w"

/*276:*/
#line 1430 "./iho.w"

void IHO::acceptSolution()
{

solution->t= trial_solution->t;

assignV(solution->y,trial_solution->y);


assignV(solution->u,trial_solution->u);

assignM(solution->S,trial_solution->S);
assignV(solution->alpha,trial_solution->alpha);

assignM(solution->A,trial_solution->A);
assignV(solution->r,trial_solution->r);

assignM(solution->Q,trial_solution->Q);
assignV(solution->rQR,trial_solution->rQR);

}





/*:276*/
#line 1680 "./iho.w"


/*236:*/
#line 327 "./iho.w"

void IHO::compTightEnclosure(interval&t_next)
{
/*238:*/
#line 348 "./iho.w"









h_trial= t_next-solution->t;



/*:238*//*240:*/
#line 404 "./iho.w"

compCoeffs();


/*:240*//*241:*/
#line 447 "./iho.w"





















interval rescale= h_trial/ad->tayl_coeff_ode->getStepsize();
hoe->getErrorTerm(predictor_excess,q+1);


scaleV(predictor_excess,v_bias::pow(rescale,q+1));

hoe->getErrorTerm(corrector_excess,p+q+1);

scaleV(corrector_excess,pow(rescale,order_trial));




/*:241*/
#line 330 "./iho.w"

/*242:*/
#line 486 "./iho.w"


/*243:*/
#line 505 "./iho.w"
















assignV(temp,solution->u);


ad->tayl_coeff_ode->set(solution->t,temp,h_trial,q);
ad->tayl_coeff_ode->compTerms();
ad->tayl_coeff_ode->sumTerms(u_next,q);






/*:243*/
#line 489 "./iho.w"


/*244:*/
#line 535 "./iho.w"

















ad->tayl_coeff_var->set(solution->t,solution->y,h_trial,q);
ad->tayl_coeff_var->compTerms();
ad->tayl_coeff_var->sumTerms(U,q);


/*:244*/
#line 492 "./iho.w"


/*245:*/
#line 561 "./iho.w"




















multMiMp(M,U,solution->A);
multMiVi(temp,M,solution->r);

multMiMp(M,U,solution->Q);
multMiVi(temp2,M,solution->rQR);

bool finite_temp= finite_interval(temp);
bool finite_temp2= finite_interval(temp2);

if(!finite_temp&&!finite_temp2)
{
control->ind= failure;
return;
}
else
{
if(finite_temp&&finite_temp2)
{
bool b= intersect(temp,temp2,temp);
assert(b);
}
else
{
if(!finite_temp)
assignV(temp,temp2);


}
}

multMiMp(M,U,solution->S);
multMiVi(temp2,M,solution->alpha);

addViVi(x,temp,temp2);





/*:245*/
#line 496 "./iho.w"


/*246:*/
#line 624 "./iho.w"






















addViVi(y_pred,u_next,predictor_excess);


addViVi(y_pred,x);


bool b= intersect(y_pred,hoe->getTrialApriori(),y_pred);
assert(b);







/*:246*/
#line 499 "./iho.w"




/*:242*/
#line 331 "./iho.w"

/*247:*/
#line 664 "./iho.w"

/*248:*/
#line 717 "./iho.w"
















setM(Fj,0.0);

for(int i= p;i>=1;i--)
{
ad->tayl_coeff_var->getTerm(M,i);

scaleM(M,C_pq[i]);
addMiMi(Fj,M);
}
addId(Fj);




/*:248*/
#line 665 "./iho.w"


/*249:*/
#line 757 "./iho.w"





















ad->tayl_coeff_var->set(t_next,y_pred,h_trial,q);

ad->tayl_coeff_var->compTerms();

setM(B,0.0);
for(int i= q;i>=1;i--)
{
ad->tayl_coeff_var->getTerm(M,i);

scaleM(M,C_qp[i]);

addMiMi(B,M);
}
addId(B);


/*:249*/
#line 668 "./iho.w"

/*250:*/
#line 795 "./iho.w"








midpoint(C,B);



/*:250*/
#line 669 "./iho.w"


/*251:*/
#line 811 "./iho.w"










bool ok= matrix_inverse->encloseMatrixInverse(Cinv,C);


if(!ok)
{
control->ind= failure;
#ifdef VNODE_DEBUG
printMessage("Could not invert the C matrix.");
#endif
return;
}


multMiMi(G,Cinv,Fj);



/*:251*/
#line 671 "./iho.w"



/*252:*/
#line 838 "./iho.w"










multMiMp(S,G,solution->S);

/*:252*/
#line 674 "./iho.w"



/*253:*/
#line 850 "./iho.w"










multMiMp(A,G,solution->A);

/*:253*/
#line 677 "./iho.w"


/*254:*/
#line 862 "./iho.w"









multMiMp(Q,G,solution->Q);









/*:254*/
#line 679 "./iho.w"


/*255:*/
#line 889 "./iho.w"













scaleV(corrector_excess,errorConstant);








/*:255*/
#line 682 "./iho.w"





/*256:*/
#line 915 "./iho.w"


/*257:*/
#line 924 "./iho.w"














setV(gj,0.0);
for(int i= p;i>=0;i--)
{
ad->tayl_coeff_ode->getTerm(term,i);
scaleV(term,C_pq[i]);
addViVi(gj,term);




}

/*:257*/
#line 917 "./iho.w"

/*258:*/
#line 951 "./iho.w"























midpoint(y_pred_point,y_pred);

assignV(temp,y_pred_point);
ad->tayl_coeff_ode->set(t_next,temp,h_trial,q);
ad->tayl_coeff_ode->compTerms();

setV(temp2,0.0);
for(int i= q;i>=0;i--)
{
ad->tayl_coeff_ode->getTerm(term,i);

scaleV(term,C_qp[i]);
addViVi(temp2,term);
}




/*:258*/
#line 918 "./iho.w"

/*259:*/
#line 993 "./iho.w"










subViVi(gj,temp2);


/*:259*/
#line 919 "./iho.w"


/*:256*/
#line 688 "./iho.w"



/*260:*/
#line 1008 "./iho.w"









addViVi(d,gj,corrector_excess);


/*:260*/
#line 691 "./iho.w"





/*261:*/
#line 1027 "./iho.w"
























subViVp(term,y_pred,y_pred_point);

multMiMi(M,Cinv,B);
subFromId(M);


multMiVi(temp,M,term);
multMiVi(w,Cinv,d);
addViVi(w,temp);



/*:261*/
#line 700 "./iho.w"


/*262:*/
#line 1064 "./iho.w"
















multMiVi(temp,A,solution->r);
multMiVi(s,Q,solution->rQR);


b= intersect(s,temp,s);
assert(b);

/*:262*/
#line 702 "./iho.w"

/*263:*/
#line 1089 "./iho.w"





























addViVi(globalExcess,s,w);
multMiVi(temp,S,trial_solution->alpha);


addViVi(temp,globalExcess);



addViVp(temp,y_pred_point);

b= intersect(trial_solution->y,y_pred,temp);
assert(b);

/*:263*/
#line 704 "./iho.w"



/*:247*/
#line 332 "./iho.w"

/*264:*/
#line 1132 "./iho.w"

trial_solution->t= t_next;





/*:264*/
#line 333 "./iho.w"

/*265:*/
#line 1141 "./iho.w"

/*266:*/
#line 1164 "./iho.w"








midpoint(trial_solution->u,trial_solution->y);

/*:266*/
#line 1142 "./iho.w"

/*267:*/
#line 1175 "./iho.w"







midpoint(trial_solution->S,S);



/*:267*/
#line 1143 "./iho.w"


/*268:*/
#line 1187 "./iho.w"



















subMiMp(S,trial_solution->S);
multMiVi(z,S,trial_solution->alpha);
addViVi(z,w);
subViVp(z,trial_solution->u);
addViVp(z,y_pred_point);


/*:268*/
#line 1145 "./iho.w"




/*269:*/
#line 1214 "./iho.w"









midpoint(trial_solution->A,A);

/*:269*/
#line 1149 "./iho.w"


/*270:*/
#line 1228 "./iho.w"


























ok= matrix_inverse->encloseMatrixInverse(Ainv,trial_solution->A);

if(ok)
{

multMiVi(temp,Ainv,z);
multMiMi(M,Ainv,A);
multMiVi(trial_solution->r,M,solution->r);
addViVi(trial_solution->r,temp);
}


/*:270*/
#line 1152 "./iho.w"


/*273:*/
#line 1335 "./iho.w"














midpoint(A_point,Q);
int n= sizeM(Q);
for(int i= 0;i<n;i++)
for(int j= 0;j<n;j++)
A_point[i][j]*= v_bias::width(solution->rQR[j]);


if(!(inf_normM(A_point)<numeric_limits<double> ::max()))
{
#ifdef VNODE_DEBUG
printMessage("The computed enclosures are too wide");
#endif
control->ind= failure;
return;
}
sortColumns(C,A_point);
b= computeQR(trial_solution->Q,C);
assert(b);


/*:273*/
#line 1154 "./iho.w"


/*271:*/
#line 1268 "./iho.w"






















b= matrix_inverse->orthogonalInverse(Ainv,trial_solution->Q);
if(b==false){
control->ind= failure;
#ifdef VNODE_DEBUG
printMessage("Could not invert the Q matrix.");
#endif
return;
}

multMiVi(temp,Ainv,z);
multMiMi(M,Ainv,Q);
multMiVi(trial_solution->rQR,M,solution->rQR);
addViVi(trial_solution->rQR,temp);


/*:271*/
#line 1157 "./iho.w"


/*272:*/
#line 1306 "./iho.w"
















multMpVi(temp,trial_solution->Q,trial_solution->rQR);
multMpVi(temp2,trial_solution->A,trial_solution->r);
if(subseteq(temp,temp2)||!ok)
{
assignM(trial_solution->A,trial_solution->Q);
assignV(trial_solution->r,trial_solution->rQR);
}




/*:272*/
#line 1159 "./iho.w"




/*:265*/
#line 334 "./iho.w"

}



/*:236*/
#line 1682 "./iho.w"

/*239:*/
#line 371 "./iho.w"


void IHO::compCoeffs()
{
if(order_trial!=control->order)
{


order_trial= control->order;
double pq= (order_trial-1)/2.0;
p= int(floor((pq)));
q= int(ceil((pq)));

assert(p+q+1==order_trial);

if(C_pq)delete[]C_pq;

C_pq= new interval[p+1];



if(C_qp)delete[]C_qp;
C_qp= new interval[q+1];


compCpq(p,q);
compCqp(p,q);
errorConstant= compErrorConstant(p,q);
}
}


/*:239*/
#line 1683 "./iho.w"

}


#line 1 "./integ.w"

/*:293*/
