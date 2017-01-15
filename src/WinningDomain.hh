/*
 * WinningDomain.hh
 *
 *  created: Jan 2017
 *   author: Matthias Rungger
 *              
 */

/** @file **/

#ifndef WINNINGDOMAIN_HH_
#define WINNINGDOMAIN_HH_

#include <vector>
#include <iostream>
#include <limits>

/* to get abs_type */
#include "TransitionFunction.hh"

#define SCOTS_WD_TYPE   "WINNINGDOMAIN"
#define SCOTS_WD_DATA   "DATA"

/** @namespace scots **/ 
namespace scots {

/**
 * @class WinningDomain
 *
 * @brief The abstract states from which the controller is winning
 *
 * The WinningDomain is computed by one of the game solving algorithms in
 * GameSovler.hh and represents the set of abstract states in {0,...,N-1} 
 * from which the controller is winning.
 *
 * For a mathematical definition see the <a href="./../../manual/manual.pdf">manual.pdf</a>.
 *
 * The class can only be \b moved \b (not copied). Also the containers to store
 * the winning states and inputs can only be \b moved into the class via the
 * constructor or the set functions. 
 *
 **/
class WinningDomain {
/* allow the write_to_file function to access the private members */
friend bool write_to_file(const WinningDomain&, const std::string&);
private:
  /** @brief size of state alphabet N **/
  abs_type m_no_states=0;
  /** @brief size of input alphabet M **/
  abs_type m_no_inputs=0;
  /** @brief used to encode loosing states **/
  abs_type m_max = std::numeric_limits<abs_type>::max();
  /**
   * @brief array of size N  \n
   * (m_winning_domain[i]=m_max if i is not winning) 
   **/
  std::vector<abs_type> m_winning_domain{};
  /**
   * @brief bool array of size N*M encoding the valid inputs\n
   *         m_inputs[i*M +j]==true iff j is a valid input at i) 
   **/
  std::vector<bool> m_inputs{}; 

public:
  /** @cond  EXCLUDE from doxygen **/
  /* default constructor */
  WinningDomain()=default;                      
  /* destructor */
  ~WinningDomain()=default;
  /* copy constructor deleted (cannot be copied) */
  WinningDomain(const WinningDomain&)=delete;
  /* move constructor */
  WinningDomain(WinningDomain&&)=default;
  /* copy assignment operator */
  WinningDomain& operator=(const WinningDomain&)=delete; 
  /* move assignment operator */
  WinningDomain& operator=(WinningDomain&&)=default;
  /* @endcond */

  /** @brief construct WinningDomain with number of states and number of abstract inputs **/
  WinningDomain(const abs_type no_states, const abs_type no_inputs) :
                m_no_states(no_states), m_no_inputs(no_inputs) {}

  /** @brief construct WinningDomain with array of winning states **/
  WinningDomain(const abs_type no_states,
                const abs_type no_inputs,
                std::vector<abs_type>&& winning_domain) :
                m_no_states(no_states), m_no_inputs(no_inputs),
                m_winning_domain(std::move(winning_domain)) {}

  /** @brief construct WinningDomain with array of winning states and valid inputs **/
  WinningDomain(const abs_type no_states,
                const abs_type no_inputs,
                std::vector<abs_type>&& winning_domain,
                std::vector<bool>&& inputs) : 
                m_no_states(no_states), m_no_inputs(no_inputs),
                m_winning_domain(std::move(winning_domain)), m_inputs(std::move(inputs)) {}
//
//  /** 
//   * @brief set the array of winning states\n
//   *        the set of states can only be set by \b moving from win_domain
//   **/
//  void set_winning_domain(std::vector<abs_type>&& win_domain) {
//		m_winning_domain = std::move(win_domain);
//  }
//
//  /** 
//   * @brief set the array of valid inputs \n
//   *        inputs can only be set by \b moving from inuts
//   **/
//  void set_inputs(std::vector<bool>&& inputs) {
//		m_inputs = std::move(inputs);
//  }
//
//  /** @brief set the size of state alphabet and input alphabet **/
//  void set_alphabet_size(const abs_type no_states, const abs_type no_inputs) {
//		m_no_states = no_states;
//		m_no_inputs = no_inputs;
//  }

  /** @brief check if state i is winning **/
  bool is_winning(const abs_type& i) {
    if(i<m_winning_domain.size() && m_winning_domain[i]!=m_max) {
      return true;
    }
    return false;
  }

  /** @brief return valid inputs associated with state i **/
  std::vector<abs_type> get_inputs(const abs_type& i) {
    std::vector<abs_type> inputs{};
    /* extract input information from m_inputs matrix */
    if(m_inputs.size()==m_no_states*m_no_inputs) {
      for(abs_type j=0; j<m_no_inputs; j++) {
        if(m_inputs[i*m_no_inputs+j]) {
          inputs.push_back(j);
        }
      }
      return inputs;
    }
    /* otherwise valid input might be written directly in m_winning_domain */
    if(is_winning(i)) {
      inputs.push_back(m_winning_domain[i]);
      return inputs;
    }
    return inputs;
  }

  /** @brief get number of winning states (=size of winning domain) **/
  abs_type get_size() {
    abs_type count=0;
    for(abs_type i=0; i<m_no_states; i++) {
      if(is_winning(i)) {
        count++;
      }
    }
    return count;
  }

  /** @brief get winning states **/
  std::vector<abs_type> get_winning_domain() {
    std::vector<abs_type> ws{};
    for(abs_type i=0; i<m_no_states; i++) {
      if(is_winning(i)) {
        ws.push_back(i);
      }
    }
    return ws;
  }

};
} /* close namespace */

#endif /* WINNINGDOMAIN_HH_ */
