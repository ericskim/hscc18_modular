/*
 * TransitionFunction.hh
 *
 *  created on: 10.01.2016
 *      author: Matthias Rungger
 * 
 */

/** 
 * @file
 * Containing the TransitionFunction and abs_type definition
 **/

#ifndef TRANSITIONFUNKTION_HH_
#define TRANSITIONFUNKTION_HH_

#include <iostream>
#include <sstream>
#include <cstring>
#include <vector>

/** @namespace scots **/ 
namespace scots {

/**
 * @brief abs_type defines type of abstract state (default = uint32_t) \n
 * determinse implicitely an upper bound on the number of states (default = 2^32-1)
 **/
using abs_type=std::uint32_t;

/**
 * @class TransitionFunction
 * 
 * @brief The transition function of the abstraction
 *
 * The transition function can either be read from file or 
 * constructed with the help of AbstractionGB.
 *  
 * The set of states is given \f$ S:=\{0,\ldots,N \}\f$ where N is stored in m_no_states
 * The set of inputs is given \f$ A:=\{0,\ldots,M \}\f$ where M is stored in m_no_inputs
 * 
 * A transition is a triple \f$ (i,j,k) \f$ where \f$i,k\in A\f$ and \f$j\in A\f$
 * - the element \f$i\f$ is referred to as pre
 * - the element \f$k\f$ is referred to as post
 * - the element \f$j\f$ is referred to as input
 * 
 * The data structure is particularly tailored to be used in the 
 * minimal and maximal fixed point computations, in which the list of pres i
 * associated with a post k and input j need to be accessed fast:\n
 * the pres itself are stored in the array m_pre\n
 * the position of the pres in m_pre is stored in the array m_pre_ptr\n
 * the number of pres is stored in  m_no_pre
 * 
 * To copy the pres of the state-input pair (k,j) to a std::vector<abs_type> pre the
 * following code is implemented in getPre(k,j)
  \verbatim
  for(std::size_t i=0; i<m_no_pre[k*M+j]; i++) {
    pre.push_back(m_pre[m_pre_ptr[k*M+j]]+i);
  }
  \endverbatim
 * 
 *
 * A transition function can only be moved and not copied.
 * 
 **/
class TransitionFunction {
public:
  /* @cond  EXCLUDE from doxygen */
  /* default constructor */
  TransitionFunction();     
  /* default destructor */
  ~TransitionFunction();    
  /* move constructor */
  TransitionFunction(TransitionFunction&&);     
  /* move asignement operator */
  TransitionFunction& operator=(TransitionFunction&&);
  /* deactivate copy constructor */
  TransitionFunction(const TransitionFunction&) = delete;     
  /* deactivate copy asignement operator */
  TransitionFunction& operator=(const TransitionFunction&) = delete;
  /* @endcond */

  /** @brief number of elements of the transition relation **/
  std::size_t getNoTransitions(void) const;  

  /** @brief return list of pre associated with action j and post k **/
  std::vector<abs_type> getPre(const abs_type& k, const abs_type& j) const; 
  /** @brief return list of post associated with action j and post i **/
  std::vector<abs_type> getPost(const abs_type& i, const abs_type& j) const;

  /** @brief number of states N **/
  abs_type m_no_states;
  /** @brief number of inputs M **/
  abs_type m_no_inputs;   
  /** @brief number of transitions T **/
  std::size_t m_no_transitions; 

  /** @brief array[N*M] containing the pre's address in the array m_pre[T] **/
  std::size_t *m_pre_ptr;    
  /** @brief array[T] containing the list of all pre */
  abs_type *m_pre;         
  /** @brief array[N*M] saving the number of post for each state-input pair (i,j) **/
  abs_type *m_no_post;      
  /** @brief array[N*M] saving the number of pre for each state-input pair (i,j) **/
  abs_type *m_no_pre;       
};

TransitionFunction::TransitionFunction() {
  m_no_states=0;
  m_no_inputs=0;
  m_no_transitions=0;

  m_pre=nullptr;
  m_pre_ptr=nullptr;
  m_no_pre=nullptr;
  m_no_post=nullptr;
}

TransitionFunction::TransitionFunction(TransitionFunction&& other) : TransitionFunction() {
  *this = std::move(other);  
}

TransitionFunction& TransitionFunction::operator=(TransitionFunction&& other) {
  m_no_states=other.m_no_states;
  m_no_inputs=other.m_no_inputs;
  m_no_transitions=other.m_no_transitions;

  m_pre=other.m_pre;
  m_pre_ptr=other.m_pre_ptr;
  m_no_pre=other.m_no_pre;
  m_no_post=other.m_no_post;

  other.m_no_states=0;
  other.m_no_inputs=0;
  other.m_no_transitions=0;

  other.m_pre=nullptr;
  other.m_pre_ptr=nullptr;
  other.m_no_pre=nullptr;
  other.m_no_post=nullptr;

  return *this;
}

TransitionFunction::~TransitionFunction() {
  delete[] m_pre;
  delete[] m_pre_ptr;
  delete[] m_no_post;
  delete[] m_no_pre;
}

std::size_t TransitionFunction::getNoTransitions(void) const {
  return m_no_transitions;
}

std::vector<abs_type> TransitionFunction::getPre(const abs_type& k, const abs_type& j) const {
  std::vector<abs_type> pre;
  if(m_pre!=nullptr) {
    if(!m_no_pre[k*m_no_inputs+j]) {
      return pre;
    }
    for(abs_type v=0; v<m_no_pre[k*m_no_inputs+j]; v++) {
      pre.push_back(m_pre[m_pre_ptr[k*m_no_inputs+j]+v]);
    }
  } else {
    std::ostringstream os;
    os << "TransitionFunction.hh: Error: Unable to get post. Transition relation is empty.";
    throw std::runtime_error(os.str().c_str());
  }
  return pre;
}

std::vector<abs_type> TransitionFunction::getPost(const abs_type& i, const abs_type& j) const {
  std::vector<abs_type> post;
  if(m_pre!=nullptr) {
    for(abs_type k=0; k<m_no_states; k++) {
      if(m_no_post[k*m_no_inputs+j]) {
        for(abs_type no=0; no<m_no_post[k*m_no_inputs+j]; no++) {
          abs_type pos=m_pre_ptr[k*m_no_inputs+j]+no;
          if(m_pre[pos]==i) {
            post.push_back(k);
          }
        }
      }
    }
  } else {
    std::ostringstream os;
    os << "TransitionFunction.hh: Error: Unable to get post. Transition relation is empty, i.e.,  m_pre and post_ are nullptr.";
    throw std::runtime_error(os.str().c_str());
  }
  return post;
}


} /* end of namespace scots */

#endif /* TRANSITIONFUNKTION_HH_ */