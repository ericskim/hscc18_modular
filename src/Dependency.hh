/*
 * Dependency.hh
 *
 *  created: Apr 2017
 *   author: Eric Kim
 */

/** @file **/
#ifndef DEPENDENCY_HH_
#define DEPENDENCY_HH_

#include <iostream>
#include <vector>
#include <unordered_set> 


#include "SymbolicSet.hh"

/** @namespace scots **/ 
namespace scots {

/**
 * @class Dependency 
 *
 * @brief Abstract base class that encodes dependencies in system dynamics
 *
 * Theoretical Background: 
 * - https://arxiv.org/abs/1704.03951
 **/
class Dependency{
protected:
  /** @brief Dimension of the state space **/
  const int state_dim;
  /** @brief Dimension of the input space **/
  const int input_dim;
  /** @brief Vector of sets of state variables located on the right hand side of update equation **/
  std::vector<std::vector<int> > state_rhs;
  /** @brief Vector of sets of state variable dependencies **/
  std::vector<std::vector<int> > state_dep;

  /** @brief Vector of sets of input variables located on the right hand side of update equation **/
  std::vector<std::vector<int> > input_rhs;
  /** @brief Vector of sets of input variable dependencies **/
  std::vector<std::vector<int> > input_dep;

public:
  Dependency(): state_dim(0), input_dim (0){}
  Dependency(int s_dim, int i_dim): state_dim(s_dim), 
                                    input_dim(i_dim), 
                                    state_rhs(s_dim, std::vector<int>()), 
                                    state_dep(s_dim, std::vector<int>()), 
                                    input_rhs(s_dim, std::vector<int>()), 
                                    input_dep(s_dim, std::vector<int>()){}

  ~Dependency() {};

  /** @brief Set right hand side of post state update equations **/
  void set_rhs(std::vector<std::vector<int> > & s_rhs, std::vector<std::vector<int> > & i_rhs){

    /** Check right hand side sizes **/
    if (s_rhs.size() != (size_t) state_dim)
      throw std::runtime_error("\nscots::Dependency state dependency vector size must be same as state dimension");
    if (i_rhs.size() !=  (size_t) state_dim)
      throw std::runtime_error("\nscots::Dependency input dependency vector size must be same as state dimension");

    /** Set state dependencies **/    
    for (size_t i= 0; i < s_rhs.size(); i++){
      if (s_rhs[i].size() == 0)
        continue;
      if ( *std::min_element(s_rhs[i].begin(), s_rhs[i].end()) < 0 || 
                 *std::max_element(s_rhs[i].begin(), s_rhs[i].end()) >= state_dim)
        throw std::runtime_error("\nscots::Dependency state index is out of bounds");
      state_rhs[i] = s_rhs[i];
    }

    /** Set input dependencies **/
    for (size_t i= 0; i < i_rhs.size(); i++){
      if (i_rhs[i].size() == 0)
        continue;
      if ( *std::min_element(i_rhs[i].begin(), i_rhs[i].end()) < 0 || 
                 *std::max_element(i_rhs[i].begin(), i_rhs[i].end()) >= input_dim)
        throw std::runtime_error("\nscots::Dependency input index is out of bounds");
      input_rhs[i] = i_rhs[i];
    }

  /** Update dependency set based off CT or DT equations **/
  set_state_dependency();
  set_input_dependency();
  }

  /** @brief Get set of state dimensions for which dimension i's post depends**/
  std::vector<int> get_state_dependency(int i){
    return state_dep[i];
  }

  /** @brief Get set of input dimensions for which dimension i's post depends**/
  std::vector<int> get_input_dependency(int i){
    return input_dep[i];
  }

  /** @brief Print out state and input dependencies **/
  friend std::ostream &operator<< (std::ostream &os, const Dependency & dep){
    for(int i = 0; i < dep.state_dim; i++){
      os << "\nPost State " << i << " Dependencies";
      os << "\nStates:";
      for (auto &j: dep.state_dep[i]){
        os << "  " << j;
      }
      os << "\nInputs:";
      for (auto &j: dep.input_dep[i]){
        os << "  " << j;
      }
    }
    os << "\n";

    return os;
  } 

  virtual void set_state_dependency() = 0; 
  virtual void set_input_dependency() = 0; 

};


/**
 * @class DT_Dependency 
 *
 * @brief Encodes dependencies with discrete time system dynamics
 **/
class DT_Dependency: public Dependency{

public:
  DT_Dependency(int s_dim, int i_dim): Dependency(s_dim, i_dim){}

  /** @brief Set DT state dependency to be same as the input right hand side**/
  void set_input_dependency(){
    for (int i = 0; i < state_dim; i++){
      input_dep[i] = input_rhs[i];
    }
  }

  /** @brief Set DT state dependency to be same as the state right hand side**/
  void set_state_dependency(){
    for (int i = 0; i < state_dim; i++){
      state_dep[i] = state_rhs[i];
    }
  } 

  ~DT_Dependency() {};

};

/**
 * @class CT_Dependency 
 *
 * @brief Encodes dependencies with continuous time system dynamics
 **/
class CT_Dependency: public Dependency{
public: 

};

} /*end namespace*/

#endif /* DEPENDENCY_HH_ */
