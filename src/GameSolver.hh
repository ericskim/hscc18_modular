/*
 * GameSolver.hh
 *
 *  created: Feb 2016
 *   author: Matthias Rungger
 * 
 */

/** @file **/

#ifndef GAMESOLVER_HH_
#define GAMESOLVER_HH_

#include <iostream>
#include <climits>
#include <stdexcept>
#include <queue>
#include <memory>
#include <utility>

#include "UniformGrid.hh"
#include "TransitionFunction.hh"
#include "WinningDomain.hh"


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
WinningDomain solve_reachability_game(const TransitionFunction& trans_function, F1& target, F2&& avoid) { 
  val_functionction_type value;
  return solve_reachability_game<F1,F2>(trans_function, target, std::forward<F2>(avoid), value);
}
/** @endcond **/

/**
 * @brief solve reachability game according to Algorithm 2 in  <a href="./../../manual/manual.pdf">manual</a>
 * 
 * @param[in] trans_function - TransitionFunction of the symbolic model
 * @param[in]  target - lambda expression of the form
 *                      \verbatim [] (abs_type &i) -> bool \endverbatim 
 *                      returns true if state i is in target set and false otherwise
 *                       
 * @param[in] avoid  - OPTIONALLY provide lambda expression of the form
 *                      \verbatim [] (abs_type &i) -> bool \endverbatim
 *                      returns true if state i is in avoid set and false otherwise
 * 
 * @param[out] value - OPTIONALLY provide val_functionction_type value to obtain the value function 
 *
 * @return -  WinningDomain that contains the set of winning states and valid inputs 
 **/
template<class F1, class F2>
WinningDomain solve_reachability_game(const TransitionFunction& trans_function, F1& target, F2&& avoid, val_functionction_type& value) {
  /* size of state alphabet */
  abs_type N=trans_function.m_no_states;
  /* size of input alphabet */
  abs_type M=trans_function.m_no_inputs;

  abs_type max = std::numeric_limits<abs_type>::max();
  if(M > max-1) {
    throw std::runtime_error("scots::solve_reachability_game: Number of inputs exceeds maximum supported value");
  }
  /* win_domain[i] = j
   * contains the input j associated with state i
   * j = max if the target is not reachable from i 
   *
   * initialize all states to loosing (win_domain[i]=max) */
  std::vector<abs_type> win_domain(N,max); 
  /* initialize value */
  value.resize(N,std::numeric_limits<double>::infinity());
  /* keep track of the number of processed post */
  abs_type* K = new abs_type[N*M];
  /* keep track of the values (corresponds to M in Alg.2)*/
  double*  edge_val = new double[N*M];

  /* init fifo */
  std::queue<abs_type> fifo;
  for(abs_type i=0; i<N; i++) {
    if(target(i) && !avoid(i)) {
      win_domain[i]=0;
      /* states in the target set are defined as loosing state */
      value[i]=SCOTS_WINNINGDOMAIN_LOSINGSTATE;
      /* states in the target are added to the fifo */
      fifo.push(i);
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
      for(abs_ptr_type v=0; v<trans_function.m_no_pre[q*M+j]; v++) {
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
          win_domain[i]=j;
        }
      }  /* end loop over all pres of state i under input j */
    }  /* end loop over all input j */
  }  /* fifo is empty */
  delete[] K;
  delete[] edge_val;

  return WinningDomain(N,M,std::move(win_domain));
}

/**
 * @brief solve invariance game according to Algorithm 1 in  <a href="./../../manual/manual.pdf">manual</a>
 * 
 * @param[in] trans_function - TransitionFunction of the symbolic model
 * @param[in] safe - lambda expression of the form
 *                    \verbatim [] (abs_type &i) -> bool \endverbatim 
 *                   returns true if state i is in safe set and false otherwise
 * @return -  WinningDomain that contains the set of winning states and valid inputs 
 **/
template<class F>
WinningDomain solve_invariance_game(const TransitionFunction& trans_function, F& safe) {
  /* size of state alphabet */
  abs_type N=trans_function.m_no_states;
  /* size of input alphabet */
  abs_type M=trans_function.m_no_inputs;
  /* helper */
//  abs_type max = std::numeric_limits<abs_type>::max(); replaced by SCOTS_WINNINGDOMAIN_LOSINGSTATE global define

  /* valid_inputs
   * boolean array of size N*M
   * valid_inputs[i*M+j]=true iff input j
   * is a valid input at state i */
  std::vector<bool> valid_inputs(N*M,false); 
  /* no_input: keep track of the number of valid inputs.
   * If no_val_in[i]==0, then state i is not winning and 
   * no_val_in[i] is set to max (see also WinningDomain) */
  std::vector<abs_type> no_val_in(N,0);
  /* keep track if an unsafe state was already added to the fifo */
  std::vector<bool> added(N,false);

  /* initialization */
  std::queue<abs_type> fifo;
  for(abs_type i=0; i<N; i++) {
    if(safe(i)) {
      for(abs_type j=0; j<M; j++) {
        if(trans_function.m_no_post[i*M+j]) {
          valid_inputs[i*M+j]=true;
          no_val_in[i]++;
        }
      }
    }
    if(!no_val_in[i]) {
      fifo.push(i);
      added[i]=true;
      /* mark no_val_in[i]=max to indicate that state i ist not winning */
      no_val_in[i]=SCOTS_WINNINGDOMAIN_LOSINGSTATE;
    }
  }

  while(!fifo.empty()) {
    abs_type k=fifo.front();
    fifo.pop();
    /* loop over all inputs */
    for(abs_type j=0; j<M; j++) {
      /* loop over all pre states of (k,j) */
      for(abs_ptr_type p=0; p<trans_function.m_no_pre[k*M+j]; p++) {
        /* (i,j,k) is a transition */
        abs_type i=trans_function.m_pre[trans_function.m_pre_ptr[k*M+j]+p];
        /* check if input j at state i is considered safe */
        if(valid_inputs[i*M+j]) {
          /* set source states i with label j as unsafe pair */
          valid_inputs[i*M+j]=false;
          no_val_in[i]--;
        }
        /* add to unsafe set if state i has no more valid inputs */
        if(!no_val_in[i] && !added[i]) {
          fifo.push(i);
          added[i]=true;
          /* mark no_val_in[i]=max to indicate that state i ist not winning */
          no_val_in[i]=SCOTS_WINNINGDOMAIN_LOSINGSTATE;
        }
      }
    }
  }
  return WinningDomain(N,M,std::move(no_val_in),std::move(valid_inputs));
}

} /* end of namespace scots */
#endif /* GAMESOLVER_HH_ */
