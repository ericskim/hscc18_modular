#ifndef STATICCONTROLLER_HH_
#define STATICCONTROLLER_HH_

#include <iostream>
#include <vector>
#include <string>
#include <limits>
#include <queue>
#include "TransitionSystem.hh"

namespace scots {

/*
 * class: StaticController
 *
 * helper class used in mexFiles to access the controller 
 * 
 *
 *
 */

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
double* val_=nullptr;

public:
/* function: StaticController
 * construction */
StaticController(TransitionSystem &ts) {
  ts_=&ts;
  N_=ts_->N_;
  M_=ts_->M_;
  domain_ = new bool[N_*M_] ();
  val_ = new double[N_] ();
}

/* function: ~StaticController
 * construction */
~StaticController() {
  delete[] domain_;
  delete[] val_;
}

/* function: getLabel
 * return the indices the controller inputs */
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

}; /* close class def */
} /* close namespace */


#endif /* STATICCONTROLLER_HH_ */
