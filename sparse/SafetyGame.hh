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


/* function: solve the safety game with respect to the set defined by
 * safe(idx)
 *
 * if safe(idx)==true -> grid point with index idx is in safe set
 * if safe(idx)==false -> grid point with index idx is not in safe set
 *
 */
//template<class F>
//void solve(F &safe) {
//        std::deque<size_t> heap;
//        /* noLabel: keep track of the number of labels leading to safe states */
//        size_t* noLabel=(size_t*)calloc(N_,sizeof(size_t));
//        /* initialization */
//        for(size_t i=0; i<N_; i++) {
//                if(safe(i)) {
//                        if(!ts_->noPost_[i]) {
//                                domain_[i]=NULL;
//                                heap.push_back(i);
//                        }
//                        else {
//                                domain_[i]=(bool*)calloc(M_,sizeof(bool));
//                                for(size_t j=0; j<M_; j++) {
//                                        if(ts_->noPost_[i][j]!=0) {
//                                                domain_[i][j]=1;
//                                                noLabel[i]++;
//                                        }
//                                }
//                        }
//                }
//                else {
//                        domain_[i]=NULL;
//                        heap.push_back(i);
//                }
//        }
//
//        while(!heap.empty()) {
//                size_t length=heap.size();
//                /* loop over all the unsafe states found in the previous loop */
//                for(size_t count=0; count<length; count++) {
//                        size_t k=heap.front();
//                        heap.pop_front();
//                        /* loop over all the labels */
//                        for(size_t j=0; j<M_; j++) {
//                                size_t numOfPre=0;
//                                if(ts_->noPre_[k])
//                                        numOfPre=ts_->noPre_[k][j];
//                                /* loop over all the source states */
//                                for(size_t v=0; v<numOfPre; v++) {
//                                        size_t pos=ts_->prePointer_[k][j]+v;
//                                        size_t i=ts_->pre_[pos];
//                                        /* check if the source states is safe */
//                                        if(noLabel[i]!=0) {
//                                                /* set source states i with label j as unsafe pair */
//                                                if(domain_[i][j]!=0) {
//                                                        domain_[i][j]=0;
//                                                        noLabel[i]--;
//                                                }
//                                                /* add to unsafe set if the source state has no label leading to safe states */
//                                                if(noLabel[i]==0) {
//                                                        heap.push_back(i);
//                                                        free(domain_[i]);
//                                                        domain_[i]=NULL;
//                                                }
//                                        }
//                                }
//                        }
//                }
//        }
//        for(size_t i=0; i<N_; i++) {
//                if(noLabel[i]!=0)
//                        val_[i]=noLabel[i];
//        }
//        free(noLabel);
//
//}     /* end of solve function*/

}; /* close class def */
} /* close namespace */

#endif /* SAFETYGAME_HH_ */
