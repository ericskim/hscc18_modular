#ifndef ABSTRACTIONGB_HH_
#define ABSTRACTIONGB_HH_

#include <iostream>
#include <cstring>
#include "UniformGrid.hh"
#include "TransitionSystem.hh"
#include "TicToc.hh"

namespace scots {
/*
 * class: AbstractionGB
 *
 * Constructs an abstraction using a growth bound according to
 * the theory in http://arxiv.org/abs/1503.03715
 *
 */
template<class state_type, class input_type>
class AbstractionGB {
private:
/* var: stateSpace_ */
const UniformGrid<state_type>* stateSpace_;
/* var: inputSpace_ */
const UniformGrid<input_type>* inputSpace_;
/* var: transitionsSystem_ */
TransitionSystem* transitionSystem_;
/* var: verbose_ */
bool verbose_;

public:
/* constructor:  AbstractionGB
 *
 * initialize with the state space <UniformGrid>  and input space  <UniformGrid>
 */
AbstractionGB(const UniformGrid<state_type> &stateSpace,
              const UniformGrid<input_type> &inputSpace,
              TransitionSystem &transitionSystem,
              bool verbose=true)
        : stateSpace_(&stateSpace),
        inputSpace_(&inputSpace),
        transitionSystem_(&transitionSystem),
        verbose_(verbose) { }

/* function:  computeTransitionRelation
 *
 * provide:
 * system_post(x,u)   -   computes the center of the hyper rectangle that
 *                        containes the attainable set
 * radius_post(r,x,u) -   computes the radius of the hyper rectangle contianing
 *                        the attainabel set
 * OPTIONAL:
 *
 * overflow(x,r)    -     check if the hyper rectangle with center x, and radius r
 *                        is an overflow symbol?
 *
 */
template<class F1, class F2>
void computeTransitionRelation(F1 &system_post, F2 &radius_post) {
  computeTransitionRelation(system_post, radius_post, [](const state_type&, const state_type&) {return false;});
}

template<class F1, class F2, class F3>
void computeTransitionRelation(F1 &system_post, F2 &radius_post, F3 &&overflow) {
  /* number of cells (=grid points) */
  size_t N=stateSpace_->getN(); 
  /* number of inputs */
  size_t M=inputSpace_->getN();
  /* state space dimension */
  int dim=stateSpace_->getDimension();
  /* number of transitions */
  size_t noT=0; 
  /* for display purpose */
  abs_type counter=0;
  /* some grid information */
  std::vector<abs_type> NN(dim);
  NN=stateSpace_->getNN();
  /* variables for managing the post */
  std::vector<abs_type> lb(dim);  /* lower-left corner */
  std::vector<abs_type> ub(dim);  /* upper-right corner */
  std::vector<abs_type> no(dim);  /* number of cells per dim */
  std::vector<abs_type> cc(dim);  /* coordinate of current cell in the post */
  /* radius of hyper interval containing the attainable set */
  state_type eta=stateSpace_->getEta();
  state_type z=stateSpace_->getZ();
  state_type r;
  /* state and input variables */
  state_type x;
  input_type u;
  /* for out of bounds check */
  state_type first;
  state_type last;
  first=stateSpace_->getFirstGridPoint();
  std::vector<abs_type> npoints(dim);
  npoints = stateSpace_->getNofGridPoints();
  for(int i=0; i<dim; i++)
    last[i]=first[i]+eta[i]*(npoints[i]-1);
  /* list of cell (pre) cell indices  */
  abs_type* pre=nullptr;
  /* number of pre indices of (i,j) */
  abs_type* noPre = new abs_type[N*M] ();  
  /* number of post indices of (i,j) */
  abs_type* noPost = new abs_type[N*M] ();  
  /* index to pre, where the cell IDs of the pre are stored */
  size_t* prePointer = new size_t[N*M];
  /* lower-left & upper-right corners of hyper rectangle of cells that cover attainable set */
  abs_type* cornerIDs = new abs_type[N*M*2];
  /* is post of (i,j) out of domain ? */
  bool* outOfDomain = new bool[N*M];
  /* loop over all cells */
  for(abs_type i=0; i<N; i++) {
    /* loop over all inputs */
    for(size_t j=0; j<M; j++) {
      outOfDomain[i*M+j]=false;
      /* cell radius (including measurement errors) */
      for(int k=0; k<dim; k++)
        r[k]=eta[k]/2.0+z[k];
      /* get center x of cell */
      stateSpace_->itox(i,x);
      /* is x an element of the overflow symbols ? */
      if(!j & overflow(x,r)) {
        for(size_t j=0; j<M; j++)
          outOfDomain[i*M+j]=true;
        break;
      }
      /* current input */
      inputSpace_->itox(j,u);
      /* integrate system and radius growth bound */
      /* the result is stored in x and r */
      radius_post(r,x,u);
      system_post(x,u);
      /* determine the cells which intersect with the attainable set: 
       * discrete hyper interval of cell indices 
       * [lb[0]; ub[0]] x .... x [lb[dim-1]; ub[dim-1]]
       * covers attainable set 
       */
      size_t npost=1;
      for(int k=0; k<dim; k++) {
        /* check for out of bounds */
        double left = x[k]-r[k]-z[k];
        double right = x[k]+r[k]+z[k];
        if(left <= first[k]-eta[k]/2.0  || right >= last[k]+eta[k]/2.0) {
          outOfDomain[i*M+j]=true;
          break;
        } 
        /* integer coordinate of lower left corner of post */
        lb[k] = static_cast<abs_type>((left-first[k]+eta[k]/2.0)/eta[k]);
        /* integer coordinate of upper right corner of post */
        ub[k] = static_cast<abs_type>((right-first[k]+eta[k]/2.0)/eta[k]);
        /* number of grid points in the post in each dimension */
        no[k]=(ub[k]-lb[k]+1);
        /* total number of post */
        npost*=no[k];
        cc[k]=0;
      }
      cornerIDs[i*(2*M)+2*j]=0;
      cornerIDs[i*(2*M)+2*j+1]=0;
      if(outOfDomain[i*M+j]) 
        continue;
      /* compute indices of post */
      for(abs_type k=0; k<npost; k++) {
        abs_type q=0;
        for(int l=0; l<dim; l++) 
          q+=(lb[l]+cc[l])*NN[l];
        cc[0]++;
        for(int l=0; l<dim-1; l++) {
          if(cc[l]==no[l]) {
            cc[l]=0;
            cc[l+1]++;
          }
        }
        /* (i,j,q) is a transition */    
        /* increment number of pres for (q,j) */ 
        noPre[q*M+j]++;
        /* store id's of lower-left and upper-right cell */
        if(k==0)
          cornerIDs[i*(2*M)+2*j]=q;
        if(k==npost-1)
          cornerIDs[i*(2*M)+2*j+1]=q;
      }
      /* increment number of transitions by number of post */
      noT+=npost;
      noPost[i*M+j]=npost;
    }
    /* print progress */
    if(verbose_ && ((double)i/(double)N*100)>counter){
      if(counter==0)
        std::cout << "1st loop: ";
      if((counter%10)==0)
        std::cout << counter;
      else if((counter%2)==0) {
        std::cout << ".";
      }
      counter++;
    }
    std::flush(std::cout); 
  }
  if(verbose_)
    std::cout << "100" << std::endl;

  /* compute prePointer */
  size_t sum=0;
  for(size_t i=0; i<N; i++) {
    for(size_t j=0; j<M; j++) {
      sum+=noPre[i*M+j];
      prePointer[i*M+j]=sum;
    }
  }
  /* allocate memory for pre list (pre[0] is used as dummy) */
  pre = new abs_type[noT];
  /* fill in pre list */
  counter=0;
  for(abs_type i=0; i<N; i++) {
    /* loop over all inputs */
    for(abs_type j=0; j<M; j++) {
    /* is x an element of the overflow symbols ? */
      if(outOfDomain[i*M+j]) 
        continue;
      /* extract lower-left and upper-bound points */
      abs_type k_lb=cornerIDs[i*2*M+2*j];
      abs_type k_ub=cornerIDs[i*2*M+2*j+1];
      abs_type npost=1;
      /* cell idx to coordinates */
      for(int k=dim-1; k>=0; k--) {
        /* integer coordinate of lower left corner */
        lb[k]=k_lb/NN[k];
        k_lb=k_lb-lb[k]*NN[k];
        /* integer coordinate of upper right corner */
        ub[k]=k_ub/NN[k];
        k_ub=k_ub-ub[k]*NN[k];
        /* number of grid points in each dimension in the post */
        no[k]=(ub[k]-lb[k]+1);
        /* totoal no of post of (i,j) */
        npost*=no[k];
        cc[k]=0;
      }
      for(abs_type k=0; k<npost; k++) {
        abs_type q=0;
        for(int l=0; l<dim; l++) 
          q+=(lb[l]+cc[l])*NN[l];
        cc[0]++;
        for(int l=0; l<dim-1; l++) {
          if(cc[l]==no[l]) {
            cc[l]=0;
            cc[l+1]++;
          }
        }
        /* (i,j,q) is a transition */
        pre[--prePointer[q*M+j]]=i;
      }
    }
    /* print progress */
    if(verbose_ && ((double)i/(double)N*100)>counter){
      if(counter==0)
        std::cout << "2nd loop: ";
      if((counter%10)==0)
        std::cout << counter;
      else if((counter%2)==0) {
        std::cout << ".";
      }
      counter++;
    }
    std::flush(std::cout); 
  }
  if(verbose_) 
    std::cout << "100" << std::endl;

  delete[] outOfDomain;
  delete[] cornerIDs;

  /* set values of abstract system */
  transitionSystem_->N_=N;
  transitionSystem_->M_=M;
  transitionSystem_->T_=noT;
  transitionSystem_->noPre_=noPre;
  transitionSystem_->noPost_=noPost;
  transitionSystem_->pre_=pre;
  transitionSystem_->prePointer_=prePointer;
}

/* function getPre */
std::vector<state_type> getPre(state_type x, input_type u) {
  abs_type k,j;
  std::vector<state_type> pre;
  std::vector<abs_type> i;
  pre.clear();
  i.clear();
  stateSpace_->xtoi(k,x);
  inputSpace_->xtoi(j,u);
  i=transitionSystem_->getPre(k,j);
  for(abs_type v=0; v<i.size(); v++) {
    state_type s;
    stateSpace_->itox(i[v],s);
    pre.push_back(s);
  }
  return pre;
}

/* function getPost */
std::vector<state_type> getPost(state_type x, input_type u) {
  abs_type i,j;
  std::vector<state_type> post;
  std::vector<abs_type> k;
  post.clear();
  k.clear();
  stateSpace_->xtoi(i,x);
  inputSpace_->xtoi(j,u);
  k=transitionSystem_->getPost(i,j);
  for(abs_type v=0; v<k.size(); v++) {
    state_type s;
    stateSpace_->itox(k[v],s);
    post.push_back(s);
  }
  return post;
}

}; /* close class def */
} /* close namespace */

#endif /* ABSTRACTIONGB_HH_ */
