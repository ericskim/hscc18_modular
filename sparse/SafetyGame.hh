#ifndef SAFETYGAME_HH_
#define SAFETYGAME_HH_

#include <iostream>
#include <limits>
#include "TransitionSystem.hh"
#include <queue>

namespace scots {
/*
 * class: SafetyGame
 *
 *
 */

/* class SafetyGame */
class SafetyGame {
friend class IO;
private:
/* var: N_
 * number of states in the transition system */
size_t N_;
/* var: M_
 * number of inputs in the transition system */
size_t M_;
/* var: ts_
 * reference to the transition system */
const TransitionSystem& ts_;
/* var: domain_
 * boolean matrix N_ x M_ 
 * domain_[i*M_+j]=true if input j 
 * is a valid input at state i */
std::unique_ptr<bool[]> domain_;

public:
/* function: SafetyGame */
SafetyGame(const TransitionSystem &ts) : ts_(ts) {
  N_=ts_.N_;
  M_=ts_.M_;
  domain_.reset(new bool[N_*M_] ());
}

/* function:  size
 * compute the number of states for which there exists a valid input */
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

/* function:  sizePairs
 * compute the number of valid state-input pairs */
abs_type sizePairs(void) const {
  abs_type n=0;
  if(!domain_)
    return n;
  for(abs_type i=0; i<N_; i++) {
    for(abs_type j=0; j<M_; j++) 
      if(domain_[i*M_+j])
        n++;
  }
  return n;
}


/* function: solve the safety game with respect to the set defined by
 * safe(idx)
 *
 * if safe(idx)==true -> grid point with index idx is in safe set
 * if safe(idx)==false -> grid point with index idx is not in safe set
 *
 */
template<class F>
void solve(F &safe) {
  /*fifo*/
  std::queue<abs_type> fifo;
  /* no_input: keep track of the number of valid inputs (that lead to safe states) */
  std::unique_ptr<abs_type[]> no_val_in(new abs_type[N_] ());
  /* keep track if an unsafe  state was already added to the fifo */
  std::unique_ptr<bool[]> added(new bool[N_] ()); 
  /* initialization */
  for(size_t i=0; i<N_; i++) {
    if(safe(i)) {
      for(size_t j=0; j<M_; j++) {
        if(ts_.noPost_[i*M_+j]) {
          domain_[i*M_+j]=1;
          no_val_in[i]++;
        }
      }
    }
    if(!no_val_in[i]) {
      fifo.push(i);
      added[i]=true;
    }
  }
  while(!fifo.empty()) {
    abs_type k=fifo.front();
    fifo.pop();
    /* loop over all inputs */
    for(abs_type j=0; j<M_; j++) {
      /* loop over all pre states of (k,j) */
      for(abs_type p=0; p<ts_.noPre_[k*M_+j]; p++) {
        /* (i,j,k) is a transition */
        abs_type i=ts_.pre_[ts_.prePointer_[k*M_+j]+p];
        /* check if input j at state i is considered safe */
        if(domain_[i*M_+j]) {
          /* set source states i with label j as unsafe pair */
          domain_[i*M_+j]=false;
          no_val_in[i]--;
        }
        /* add to unsafe set if state i has no more valid inputs */
        if(!no_val_in[i] && !added[i]) {
          fifo.push(i);
          added[i]=true;
        }
      }
    }
  }
}/* end of solve function*/

}; /* close class def */
} /* close namespace */

#endif /* SAFETYGAME_HH_ */
