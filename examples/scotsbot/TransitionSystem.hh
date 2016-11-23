/*
 * TransitionSystem.hh
 *
 *  created on: 10.01.2016
 *      author: rungger
 */

#ifndef TRANSITIONSYSTEM_HH_
#define TRANSITIONSYSTEM_HH_

#include <iostream>
#include <cstring>


namespace scots {

/* forward declaration of abstraction growth bound class */
template<class stateType, class inputType> class AbstractionGB;
typedef unsigned int state_size_t;
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
  template<class stateType, class inputType>
  friend class AbstractionGB;
  friend class ReachabilityGame;
  friend class SafetyGame;
  friend class IO;
  friend class Game;
  private:
    /* var: N_
     * number of states */
    size_t N_;
    /* var: M_
     * number of inputs */
    size_t M_;
    /* var: T_
     * number of transitions */
    size_t T_;

    /* var: pre_
     * save the pre states of each state-action pair */
    state_size_t *pre_;
    /* var: post_
     * save the post states of each state-action pair */
    state_size_t *post_;
    /* var: postPointer_
     * 2-dimensional array saving the addresses for post */
    state_size_t **postPointer_;
    /* var: prePointer_
     * 2-dimensional array saving the addresses for pre */
    state_size_t **prePointer_;
    /* var: noPost_
     * 2-dimensional array saving the number of post for each state-action pair */
    short **noPost_;
    /* var: noPre_
     * 2-dimensional array saving the number of pre for each state-action pair */
    short **noPre_;

  public:
    /* constructor: TransitionSystem
     *
     * initialize all the pointer and variables to zero
     *
     */
    TransitionSystem() {
      N_=0;
      M_=0;
      pre_=0;
      post_=0;
      T_=0; /* number of transitions */
      prePointer_=NULL;
      postPointer_=NULL;
      noPost_=NULL;
      noPre_=NULL;
//      avgNoPost_=0; /* avarage number of elements in the post F(x,u)*/
      
    }
	
    ~TransitionSystem(){
      free(pre_);
      pre_=NULL;
      free(post_);
      post_=NULL;
      for(size_t i=0;i<N_;i++) {
        if(prePointer_[i])
          free(prePointer_[i]);
        if(postPointer_[i])
          free(postPointer_[i]);
        if(noPre_[i])
          free(noPre_[i]);
        if(noPost_[i])
          free(noPost_[i]);
      }
      free(prePointer_);
      free(postPointer_);
      free(noPre_);
      free(noPost_);

    }

    size_t getNoTransitions(void) {
      return T_;
    }

    /* function: getPre */
    std::vector<size_t> getPre(size_t k, size_t j) {
      std::vector<size_t> vec;
      vec.clear();
      if(pre_!=NULL) {
        if(!noPre_[k])
          return vec;
        size_t num=noPre_[k][j];
        size_t pos=0;
        for(size_t v=0;v<num;v++) {
          pos=prePointer_[k][j]+v;
          vec.push_back(pre_[pos]);
        }
      }
      else if(post_!=NULL){
        for(size_t i=0;i<N_;i++) {
          if(noPost_[i]) {
            if(noPost_[i][j]!=0 && postPointer_[i]!=NULL) {
              for(short no=0;no<noPost_[i][j];no++) {
                size_t pos=postPointer_[i][j]+no;
                if(post_[pos]==k)
                  vec.push_back(i);
              }
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
    std::vector<size_t> getPost(size_t i, size_t j) {
      std::vector<size_t> vec;
      vec.clear();
      if(post_!=NULL) {
        if(!post_[i])
          return vec;
        size_t num=noPost_[i][j];
        size_t pos=0;
        for(size_t v=0;v<num;v++) {
          pos=postPointer_[i][j]+v;
          vec.push_back(post_[pos]);
        }
      }
      else if(pre_!=NULL) {
        for(size_t k=0;k<N_;k++) {
          if(noPre_[k]) {
            if(noPre_[k][j]!=0) {
              for(short no=0;no<noPre_[k][j];no++) {
                size_t pos=prePointer_[k][j]+no;
                if(pre_[pos]==i)
                  vec.push_back(k);
              }
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

//    size_t getNoAvgPost(void) {
//      return avgNoPost_;
//    }
}; /* close class def */
} /* close namespace */

#endif /* TRANSITIONSYSTEM_HH_ */
