/*226:*/
#line 832 "./hoe.w"

#ifndef HOE_H
#define HOE_H
namespace vnodelp{
/*204:*/
#line 83 "./hoe.w"

class HOE
{
public:
HOE(int n);

void compAprioriEnclosure(const
interval&t0,const iVector&y0,bool&info);
void acceptSolution();
/*221:*/
#line 729 "./hoe.w"

void set(Control*ctrl,AD*ad)
{
control= ctrl;
tayl_coeff= ad->tayl_coeff_ode;
}

void init(const interval&t0,const iVector&y0){apriori->init(t0,y0);}
void setTrialStepsize(double h0)
{
h_trial= h0;
}

void setTrialOrder(int order0)
{
order_trial= order0;
}




/*:221*/
#line 92 "./hoe.w"

/*222:*/
#line 754 "./hoe.w"



double getStepsize()const{return h;}

double getTrialStepsize()const
{
return h_trial;
}

const interval&getT()const
{
return apriori->t;
}

interval getTrialT()const
{
return apriori_trial->t;
}


const iVector&getApriori()const
{
return apriori->y;
}

const iVector&getTrialApriori()const
{
return apriori_trial->y;
}



/*:222*//*223:*/
#line 791 "./hoe.w"






void getErrorTerm(iVector&e,int i)const
{
tayl_coeff->getTerm(e,i);
}



/*:223*/
#line 93 "./hoe.w"


~HOE();




private:
Apriori*apriori_trial,*apriori;

Control*control;
AD_ODE*tayl_coeff;

double h,h_next,h_trial,t_trial;
int order_trial;


iVector term,p,u,v;
const interval one;




interval comp_beta(const iVector&v,const iVector&u,int k);
};



/*:204*/
#line 836 "./hoe.w"

}
#endif


/*:226*/
