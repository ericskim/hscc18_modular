/*
 * TransitionSystem.hh
 *
 *  created on: 10.01.2016
 *      author: rungger
 */

// short changed to unsigned short

#ifndef TRANSITIONSYSTEM_HH_
#define TRANSITIONSYSTEM_HH_

#include <iostream>
#include <cstring>
#include <cstdint>


namespace scots {

/* forward declaration of abstraction growth bound class */
template<class state_type, class input_type> class AbstractionGB;

/* define type of abstract state
 * determines an upper limit on the number of states
  */
using abs_type=uint32_t;

/*
 * class: TransitionSystem
 *
 * Contains the transition relation of a transition system.
 * The data structure is particularly tailored to be used in the
 * minimal and maximal fixed point computations.
 *
 */
class TransitionSystem {
/* friend classes */
friend class IO;
template<class state_type, class input_type>
friend class AbstractionGB;
friend class StaticController;
friend class ReachabilityGame;
friend class SafetyGame;
private:
/* var: N_
 * number of states */
abs_type N_=0;
/* var: M_
 * number of inputs */
abs_type M_=0;
/* var: T_
 * number of transitions */
size_t T_=0;

/* var: pre_
 * list of pre cell ids */
abs_type *pre_=nullptr;
/* var: post_
 * list of post cell ids */
abs_type *post_=nullptr;
/* var: postPointer_
 * array saving the addresses for post */
size_t *postPointer_=nullptr;
/* var: prePointer_
 * array saving the addresses for pre */
size_t *prePointer_=nullptr;
/* var: noPost_
 * array saving the number of post for each state-action pair */
abs_type *noPost_=nullptr;
/* var: noPre_
 * array saving the number of pre for each state-action pair */
abs_type *noPre_=nullptr;

public:
~TransitionSystem(){
  if(pre_)
    delete[] pre_;
  if(post_)
    delete[] pre_;
  if(prePointer_)
    delete[] prePointer_;
  if(postPointer_)
    delete[] postPointer_;
  if(noPre_)
    delete[] noPre_;
  if(noPost_)
    delete[] noPost_;
}

size_t getNoTransitions(void) {
  return T_;
}

/* function: getPre */
std::vector<abs_type> getPre(abs_type k, abs_type j) {
  std::vector<abs_type> vec;
  vec.clear();
  if(pre_!=nullptr) {
    if(!noPre_[k*M_+j])
      return vec;
    abs_type num=noPre_[k*M_+j];
    abs_type pos=0;
    for(abs_type v=0; v<num; v++) {
      pos=prePointer_[k*M_+j]+v;
      vec.push_back(pre_[pos]);
    }
  }
  else if(post_!=nullptr) {
    for(abs_type i=0; i<N_; i++) {
      if(noPost_[i*M_+j]) {
        for(abs_type no=0; no<noPost_[i*M_+j]; no++) {
          abs_type pos=postPointer_[i*M_+j]+no;
          if(post_[pos]==k)
            vec.push_back(i);
        }
      }
    }
  }
  else {
    std::ostringstream os;
    os << "TransitionSystem.hh: Error: Unable to get post. Transition relation is empty, i.e.,  pre_ and post_ are NULL.";
    throw std::runtime_error(os.str().c_str());
  }
  return vec;
}

/* function: getPost */
std::vector<abs_type> getPost(abs_type i, abs_type j) {
  std::vector<abs_type> vec;
  vec.clear();
  if(post_!=nullptr) {
    if(!noPost_[i*M_+j])
      return vec;
    abs_type num=noPost_[i*M_+j];
    abs_type pos=0;
    for(abs_type v=0; v<num; v++) {
      pos=postPointer_[i*M_+j]+v;
      vec.push_back(post_[pos]);
    }
  }
  else if(pre_!=nullptr) {
    for(abs_type k=0; k<N_; k++) {
      if(noPre_[k*M_+j]) {
        for(abs_type no=0; no<noPre_[k*M_+j]; no++) {
          abs_type pos=prePointer_[k*M_+j]+no;
          if(pre_[pos]==i)
            vec.push_back(k);
        }
      }
    }
  }
  else {
          std::ostringstream os;
          os << "TransitionSystem.hh: Error: Unable to get post. Transition relation is empty, i.e.,  pre_ and post_ are NULL.";
          throw std::runtime_error(os.str().c_str());
  }
  return vec;
}
}; /* close class def */
} /* close namespace */

#endif /* TRANSITIONSYSTEM_HH_ */
