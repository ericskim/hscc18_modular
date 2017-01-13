/*
 * GameSolver.hh
 *
 *  created on: 12.02.2016
 *      author: Matthias Rungger
 * 
 */

/** @file GameSolver.hh **/

#ifndef GAMESOLVER_HH_
#define GAMESOLVER_HH_

#include <iostream>
#include <climits>
#include <stdexcept>
#include <queue>
#include <memory>

#include "TransitionFunction.hh"
#include "WinningDomain.hh"
//#include "FileHandler.hh"



/** @namespace scots **/ 
namespace scots {

/**
 * @brief val_functionction_type (value function type) is used in solve_reachability_game \n
 * to represent the value function 
 **/
using val_functionction_type = std::vector<double>;

/** @cond **/
template<class F>
WinningDomain solve_reachability_game(const TransitionFunction& trans_function, F& target) { 
  val_functionction_type value;
  return solve_reachability_game<F>(trans_function, target, [](const abs_type&) noexcept {return false;}, value);
}

template<class F>
WinningDomain solve_reachability_game(const TransitionFunction& trans_function, F& target, val_functionction_type& value) { 
  return solve_reachability_game<F>(trans_function, target, [](const abs_type&) noexcept {return false;}, value);
}

template<class F1, class F2>
WinningDomain solve_reachability_game(const TransitionFunction& trans_function, F1& target, F2& avoid) { 
  val_functionction_type value;
  return solve_reachability_game<F1,F2>(trans_function, target, avoid, value);
}
/** @endcond **/

/**
 * @brief solve reachability game according to Algorithm 2 in  <a href="./../../manual/manual.pdf">manual</a>
 * 
 * @param[in] trans_function - TransitionFunction of the symbolic model
 * @param[in]  target - lambda function with signature  
 *                      \verbatim [] (abs_type &i) -> bool \endverbatim 
 *                      returns true if state i is in target set and false otherwise
 *                       
 * @param[in] avoid  - OPTIONALLY provide lambda function with signature 
 *                      \verbatim [] (abs_type &i) -> bool \endverbatim
 *                      returns true if state i is in avoid set and false otherwise
 * 
 * @param[out] value - OPTIONALLY provide val_functionction_type value to obtain the value function 
 *
 * @return -  WinningDomain that contains the set of winning states and valid inputs 
 *
 **/
template<class F1, class F2>
WinningDomain solve_reachability_game(const TransitionFunction& trans_function, F1& target, F2&& avoid, val_functionction_type& value) {
  /* size of state alphabet */
  std::size_t N=trans_function.m_no_states;
  /* size of input alphabet */
  std::size_t M=trans_function.m_no_inputs;
  /* 
   * winning_states[i] = j
   * contains the input j associated with state i
   * j = max if the target is not reachable from i 
   */
  abs_type max = std::numeric_limits<abs_type>::max();
  std::cout << "max max " << max << std::endl;
  if(M > max-1) {
    throw std::runtime_error("scots::solve_reachability_game: Number of inputs exceeds maximum supported value");
  }
  /* initialize all states to loosing (winning_states[i]=max) */
  std::vector<abs_type> winning_states(N,max); 

  /* initialize value */
  value.resize(N,std::numeric_limits<double>::infinity());
  /* use fifo list */
  std::queue<abs_type> fifo;
  /* keep track of the number of processed post */
  std::unique_ptr<abs_type[]> K (new abs_type[N*M]);
  /* keep track of the values (corresponds to M in Alg.2)*/
  std::unique_ptr<double[]> edge_val (new double[N*M]);
  /* keep track of number of valid inputs */
  abs_type counter = 0;

  /* init fifo */
  for(abs_type i=0; i<N; i++) {
    if(target(i) && !avoid(i)) {
      winning_states[i]=0;
      /* states in the target set have value zero */
      value[i]=0; 
      /* states in the target are added to the fifo */
      fifo.push(i);
      counter++;
    }
    for(abs_type j=0; j<M; j++) {
      edge_val[i*M+j]=0;
      K[i*M+j]=trans_function.m_no_post[i*M+j];
    }
  }

  /* main loop */
  while(!fifo.empty()) {
    /* get state to be processed */
    abs_type q=fifo.front();
    fifo.pop();
    /* loop over each input */
    for(abs_type j=0; j<M; j++) {
      /* loop over pre's associated with this input */
      for(abs_type v=0; v<trans_function.m_no_pre[q*M+j]; v++) {
        abs_type i=trans_function.m_pre[trans_function.m_pre_ptr[q*M+j]+v];
        if(avoid(i))
          continue;
        /* (i,j,q) is a transition */
        /* update the number of processed posts */
        K[i*M+j]--;
        /* update the max value of processed posts */
        edge_val[i*M+j]=(edge_val[i*M+j]>=1+value[q] ? edge_val[i*M+j] : 1+value[q]);
        /* check if for node i and input j all posts are processed */
        if(!K[i*M+j] && value[i]>edge_val[i*M+j]) {
          fifo.push(i);
          value[i]=edge_val[i*M+j]; 
          winning_states[i]=j;
          counter++;
        }
      }  /* end loop over all pres of state i under input j */
    }  /* end loop over all input j */
  }  /* fifo is empty */

  return WinningDomain(N,M,std::move(winning_states));
}

//template<class F>
//WinningDomain solve_invariance_game(const TransitionFunction& trans_function, F& safe) {
//  /* var: N
//   * number of states in the transition system */
//  size_t N=ts_.m_NrOfStates;
//  /* var: M
//   * number of inputs in the transition system */
//  size_t M=ts_.m_NrOfInputs;
//  /* var: domain_
//   * boolean matrix N x M
//   * domain_[i*M+j]=true if input j
//   * is a valid input at state i */
//  nonDeterministicMap domain_; ///
//  domain_.resize(N,M); ///
//
//  std::queue<abs_type> fifo;
//  /* no_input: keep track of the number of valid inputs (that lead to safe states) */
//  std::unique_ptr<abs_type[]> no_val_in(new abs_type[N] ());
//  /* keep track if an unsafe  state was already added to the fifo */
//  std::unique_ptr<bool[]> added(new bool[N] ());
//
//  /* initialization */
//  for(size_t i=0; i<N; i++) {
//    if(safe(i)) {
//      for(size_t j=0; j<M; j++) {
//        if(ts_.noPost_[i*M+j]) {
//          domain_.set(i,j,true); ///
//          no_val_in[i]++;
//        }
//      }
//    }
//    if(!no_val_in[i]) {
//      fifo.push(i);
//      added[i]=true;
//    }
//  }
//  while(!fifo.empty()) {
//    abs_type k=fifo.front();
//    fifo.pop();
//    /* loop over all inputs */
//    for(abs_type j=0; j<M; j++) {
//      /* loop over all pre states of (k,j) */
//      for(abs_type p=0; p<ts_.noPre_[k*M+j]; p++) {
//        /* (i,j,k) is a transition */
//        abs_type i=ts_.pre_[ts_.prePointer_[k*M+j]+p];
//        /* check if input j at state i is considered safe */
//        if(domain_.get(i,j))
//        {
//          /* set source states i with label j as unsafe pair */
//          domain_.set(i,j,false);
//          no_val_in[i]--;
//        }
//        /* add to unsafe set if state i has no more valid inputs */
//        if(!no_val_in[i] && !added[i]) {
//          fifo.push(i);
//          added[i]=true;
//        }
//      }
//    }
//  }
//  if(minimalPairs) ///
//  {
//    *minimalPairs = domain_.extractMinimalValidPairs();
//  }
//
//  return domain_;
//}

} /* end of namespace scots */

#endif /* GAMESOLVER_HH_ */
