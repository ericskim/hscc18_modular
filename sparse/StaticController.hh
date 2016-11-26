#ifndef STATICCONTROLLER_HH_
#define STATICCONTROLLER_HH_

#include <iostream>
#include <vector>
#include <string>
#include <limits>
#include <queue>
#include "TransitionSystem.hh"

namespace scots {

class StaticController {
friend class IO;
/* var: N_
 * number of states in the transition system */
size_t N_;
/* var: M_
 * number of labels in the transition system */
size_t M_;
/* var: ts_
 * pointer to the transition system */
 TransitionSystem *ts_=nullptr;
/* var: domain_
 * contains the controller domain N_ x M_ */
bool* domain_=nullptr;
/* var: val_ 
 * contains the value function */
float* val_=nullptr;

public:
/* function: StaticController
 * construction */
StaticController( TransitionSystem &ts) {
  ts_=&ts;
  N_=ts_->N_;
  M_=ts_->M_;
  domain_ = new bool[N_*M_] ();
  val_ = new float[N_] ();
}

/* function: ~StaticController
 * construction */
~StaticController() {
  delete[] domain_;
  delete[] val_;
}

/* function: getLabel
 * return the index of the controller input */
std::vector<abs_type> getLabel(abs_type idx) const {
  std::vector<abs_type> label;
  if(!domain_)
    return label;
  for (abs_type j=0; j<M_; j++) {
    if (domain_[idx*M_+j])
      label.push_back(j);
    return label;
  }
  return label;
}

/* function:  size
 * compute the number of states for which there exists a valid label/input value */
abs_type size(void) const {
  abs_type n=0;
  if(!domain_)
    return n;
  for(abs_type i=0; i<N_; i++) {
    for(abs_type j=0; j<M_; j++) {
      if(domain_[i*M_+j]){
        n++;
        break;
      }
    }
  }
  return n;
}

/* function: getStateIDs
 * copy the domain indices to an integer vector of length
 * numberOfStatesInDomain() */
void getStateIDs(abs_type *idx) const {
  if(!domain_)
    return;
  for(abs_type j=0, i=0; i<N_; i++) {
    for(abs_type k=0; k<M_; k++) {
      if(domain_[i*M_+k]) {  
        idx[j]=i;
        break;
        j++;
      }
    }
  }
}

/* function: getDomainIdx
 * return the domain indices to an vector */
std::vector<abs_type> getDomainIdx() const {
  std::vector<abs_type> idx;
  if(!domain_)
    return idx;
  for(abs_type i=0; i<N_; i++) {
    for(abs_type j=0; j<M_; j++) {
      if(domain_[i*M_+j]) {  
        idx.push_back(i);
        break;
      }
    }
  }
  return idx;
}

/* function: getValue
 * copy the values of the value function to a double vector of length
 * sizeOfDomain() */
void getValue(double *val) const {
  if(!domain_ || !val_)
    return;
  for(abs_type j=0, i=0; i<N_; i++) {
    for(abs_type k=0; k<M_; k++) {
      if(domain_[i*M_+k]) {  
        val[j]=val_[i];
        break;
        j++;
      }
    }
  }
}

/* function: inDomain
 * does there exist a valid input for the state with index idx */
bool inDomain(abs_type idx) const {
  for(abs_type j=0; j<M_; j++) {
    if(domain_[idx*M_+j])
      return true;
  }
  return false;
}

/* function: reach 
 * solve the reachability problem with respect to target set
 *
 * if target(idx)==true -> grid point with index idx is in target set
 * if target(idx)==false -> grid point with index idx is not in target set
 *
 */
template<class F>
void reach(F &target) {
  /* use fifo list */
  std::queue<abs_type> fifo;
  /* controller */
  abs_type* label = new abs_type[N_];
  /* keep track of the number of processed post */
  //abs_type* k = new abs_type[N_*M_];
  /* keep track of the values */
  float* edge_val = new float[N_*M_];

  /* init fifo */
  for(abs_type i=0; i<N_; i++) {
    val_[i]=std::numeric_limits<float>::infinity();
    if(target(i)) {
      domain_[i*M_]=true;
      /* nodes in the target set have value zero */
      val_[i]=0;
      //fifo[last++]=i;
      fifo.push(i);
    }
    for(abs_type j=0; j<M_; j++) 
      edge_val[i*M_+j]=0;
  }

  while(!fifo.empty()) {
    /* get state to be processed */
    abs_type q=fifo.front();
    fifo.pop();
    /* save input label to domain_ */
    domain_[q*M_+label[q]]=true;
    /* loop over each label */
    for(abs_type j=0; j<M_; j++) {
      /* loop over pre's associated with this label */
      for(abs_type v=0; v<ts_->noPre_[q*M_+j]; v++) {
        abs_type i=ts_->pre_[ts_->prePointer_[q*M_+j]+v];
        /* (i,j,q) is a transition */
        /* update the number of processed posts */
        ts_->noPost_[i*M_+j]--;
        /* update the max value of processed posts */
        edge_val[i*M_+j]=(edge_val[i*M_+j]>=1+val_[q] ? edge_val[i*M_+j] : 1+val_[q]);
        /* check if for node i and label j all posts are processed */
        if(!ts_->noPost_[i*M_+j] && val_[i]>edge_val[i*M_+j]) {
          fifo.push(i);
          val_[i]=edge_val[i*M_+j];
          label[i]=j;
        }  
      }  /* end loop over all pres of node i  under label j */
    }  /* end loop over all label j */
  }  /* fifo is empty */

  delete[] label;
  delete[] edge_val;
}

}; /* close class def */
} /* close namespace */


#endif /* STATICCONTROLLER_HH_ */
