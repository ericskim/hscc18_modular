/*
 * ReachabilityGame.hh
 *
 *  created on: 08.01.2016
 *      author: rungger
 */
#ifndef REACHABILITYGAME_HH_
#define REACHABILITYGAME_HH_

#include <iostream>
#include <climits>
#include <stdexcept>
#include <queue>

#include "TransitionSystem.hh"


namespace scots {

/*
 * class: ReachabilityGame
 *
 *
 * Implementation of the Dijkstra algorithm for hyper-graphs according to
 *
 * Gallo G, Longo G, Pallottino S, Nguyen S. Directed hypergraphs and
 * applications. Discrete applied mathematics. 1993 Apr 27;42(2):177-201.
 *
 * to solve a reachability problem
 *
 */
class ReachabilityGame {
friend class IO;
/* var: N_
 * number of states in the transition system */
size_t N_;
/* var: M_
 * number of inputs in the transition system */
size_t M_;
/* var: ts_
 * pointer to the transition system */
TransitionSystem *ts_=nullptr;
/* var: inputs_
 * inputs_[i] = j 
 * contains the input j associated with state i
 * j=-1 if the target is not reachable from i */
int* inputs_=nullptr;
/* var: val_ 
 * contains the value function */
double* val_=nullptr;

public:
/* function: ReachabilityGame
 * construction */
ReachabilityGame( TransitionSystem &ts) {
  ts_=&ts;
  N_=ts_->N_;
  M_=ts_->M_;
  if(M_> INT_MAX)
    throw std::runtime_error("scots::ReachabilityGame: Number of inputs exceeds maximum supported value");
  inputs_ = new int[N_];
  val_ = new double[N_];
  for(abs_type i=0; i<N_; i++)
    inputs_[i]=-1;
}

/* function: ~ReachabilityGame
 * construction */
~ReachabilityGame() {
  delete[] inputs_;
  delete[] val_;
}

/* function:  size
 * compute the number of states for which there exists a valid input value */
abs_type size(void) const {
  abs_type n=0;
  if(!inputs_)
    return n;
  for(abs_type i=0; i<N_; i++) {
    if(inputs_[i]!=-1)
      n++;
  }
  return n;
}

/* function: solve 
 * solve the reachability problem with respect to target set
 *
 * if target(idx)==true -> grid point with index idx is in target set
 * if target(idx)==false -> grid point with index idx is not in target set
 *
 */
template<class F>
void solve(F &target) {
  /* use fifo list */
  std::queue<abs_type> fifo;
  /* keep track of the number of processed post */
  abs_type* K = new abs_type[N_*M_];
  /* keep track of the values */
  float* edge_val = new float[N_*M_];

  /* init fifo */
  for(abs_type i=0; i<N_; i++) {
    val_[i]=std::numeric_limits<float>::infinity();
    if(target(i)) {
      inputs_[i]=0;
      /* nodes in the target set have value zero */
      val_[i]=0;
      fifo.push(i);
    }
    for(abs_type j=0; j<M_; j++) {
      edge_val[i*M_+j]=0;
      K[i*M_+j]=ts_->noPost_[i*M_+j];
    }
  }

  while(!fifo.empty()) {
    /* get state to be processed */
    abs_type q=fifo.front();
    fifo.pop();
    /* loop over each input */
    for(abs_type j=0; j<M_; j++) {
      /* loop over pre's associated with this input */
      for(abs_type v=0; v<ts_->noPre_[q*M_+j]; v++) {
        abs_type i=ts_->pre_[ts_->prePointer_[q*M_+j]+v];
        /* (i,j,q) is a transition */
        /* update the number of processed posts */
        K[i*M_+j]--;
        /* update the max value of processed posts */
        edge_val[i*M_+j]=(edge_val[i*M_+j]>=1+val_[q] ? edge_val[i*M_+j] : 1+val_[q]);
        /* check if for node i and input j all posts are processed */
        if(!K[i*M_+j] && val_[i]>edge_val[i*M_+j]) {
          fifo.push(i);
          val_[i]=edge_val[i*M_+j];
          inputs_[i]=j;
        }  
      }  /* end loop over all pres of state i under input j */
    }  /* end loop over all input j */
  }  /* fifo is empty */


  delete[] edge_val;
  delete[] K;
}

}; /* close class def */
} /* close namespace */

#endif
