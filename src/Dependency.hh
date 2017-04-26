/*
 * Dependency.hh
 *
 *  created: Apr 2017
 *   author: Eric Kim
 */

/** @file **/

#include <iostream>
#include <vector>
#include <unordered_set> 


#include "SymbolicSet.hh"


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
  Dependency(int s_dim, int i_dim): state_dim(s_dim), input_dim(i_dim), state_rhs(s_dim, {}), state_dep(s_dim, {}), input_rhs(i_dim, {}), input_dep(i_dim, {}){
  }
  ~Dependency();

  /** @brief Set right hand side of equations **/
  void set_rhs(std::vector<std::vector<int> > & s_rhs, std::vector<std::vector<int> > & i_rhs){
    /** Set state dependencies **/
    if (s_rhs.size() != (size_t) state_dim)
      throw std::runtime_error("\nscots::Dependency RHS vector size does not match state dimension");
    for (size_t i= 0; i < s_rhs.size(); i++){
      if ( *std::min_element(s_rhs[i].begin(), s_rhs[i].end()) < 0 || 
                 *std::max_element(s_rhs[i].begin(), s_rhs[i].end()) >= state_dim)
        throw std::runtime_error("\nscots::Dependency index is out of bounds");
      state_rhs[i] = s_rhs[i];
    }

    /** Set input dependencies **/
    if (i_rhs.size() !=  (size_t) input_dim)
      throw std::runtime_error("\nscots::Dependency RHS vector size does not match input dimension");
    for (size_t i= 0; i < i_rhs.size(); i++){
      if ( *std::min_element(i_rhs[i].begin(), i_rhs[i].end()) < 0 || 
                 *std::max_element(i_rhs[i].begin(), i_rhs[i].end()) >= state_dim)
        throw std::runtime_error("\nscots::Dependency index is out of bounds");
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

  virtual void set_state_dependency() = 0; 
  virtual void set_input_dependency() = 0; 

};


/**
 * @class DT_Dependency 
 *
 * @brief Encodes dependencies with discrete time system dynamics
 **/
class DT_Dependency: public Dependency{

  DT_Dependency(int s_dim, int i_dim): Dependency(s_dim, i_dim){}
  void set_state_dependency(){
    for (int i = 0; i < state_dim; i++){
      state_dep[i] = state_rhs[i];
    }
  } 
  void set_input_dependency(){
    for (int i = 0; i < state_dim; i++){
      input_dep[i] = input_rhs[i];
    }
  }
  public:
};

/**
 * @class DT_Dependency 
 *
 * @brief Encodes dependencies with continuous time system dynamics
 **/
class CT_Dependency: public Dependency{
public: 

};

}
