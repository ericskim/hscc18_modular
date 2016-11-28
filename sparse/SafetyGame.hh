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
 * pointer to the transition system */
TransitionSystem *ts_=nullptr;
/* var: domain_
 * boolean matrix N_ x M_ 
 * domain_[i*M_+j]=true if input j 
 * is a valid input at state i */
bool* domain_=nullptr;

public:
/* function: SafetyGame */
SafetyGame(TransitionSystem &ts) {
  ts_=&ts;
  N_=ts_->N_;
  M_=ts_->M_;
  domain_ = new bool[N_*M_] ();
}

/* function: ~SafetyGame */
~SafetyGame() {
  delete[] domain_;
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
  std::queue<abs_type> fifo;

  /* no_input: keep track of the number of valid inputs (that lead to safe states) */
  abs_type* no_val_in= new abs_type[N_] ();
  /* keep track if an unsafe  state was already added to the fifo */
  bool* added = new bool[N_] ();

  /* initialization */
  for(size_t i=0; i<N_; i++) {
    if(safe(i)) {
      for(size_t j=0; j<M_; j++) {
        if(ts_->noPost_[i*M_+j]) {
          domain_[i*M_+j]=1;
          no_val_in[i]++;
        }
      }
    }
    if(!no_val_in[i] && !added[i]) {
      fifo.push(i);
      added[i]=true;
    }
  }

  while(!fifo.empty()) {
    abs_type length=fifo.size();
    /* loop over all the unsafe states found in the previous loop */
    for(abs_type count=0; count<length; count++) {
      abs_type k=fifo.front();
      fifo.pop();
      /* loop over all inputs */
      for(abs_type j=0; j<M_; j++) {
        /* loop over all pre states of (k,j) */
        for(abs_type p=0; p<ts_->noPre_[k*M_+j]; p++) {
          /* (i,j,k) is a transition */
          abs_type i=ts_->pre_[ts_->prePointer_[k*M_+j]+p];
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
  }
  delete[] no_val_in;
  delete[] added;
}/* end of solve function*/

}; /* close class def */
} /* close namespace */

#endif /* SAFETYGAME_HH_ */
