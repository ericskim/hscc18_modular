/*
 *  FixedPoint.hh
 *
 *  created: Jan 2017
 *   author: Matthias Rungger
 */

/**
 * @file
 * @brief contains the EnfPre class and some fixed point algorithms for synthesis 
 **/

#ifndef FIXEDPOINT_HH_
#define FIXEDPOINT_HH_

#include <iostream>
#include <memory>

#include "SymbolicSet.hh"
#include "SymbolicModel.hh"


namespace scots {

/**
 * @class EnfPre
 * 
 * @brief computes of the enforcable predecessor 
 * 
 * Let \f$ F:X\times U \rightrightarrows X\f$ be the transition function and \f$
 * Z\subseteq X\f$, then \n\n
 *
 *  \f$\mathrm{pre}(Z)=\{ (x,u) \in X\times U \mid F(x,u)\neq \emptyset \wedge  F(x,u)\subseteq Z\}\f$
 *
 **/
class EnfPre {
private:
  /* stores the permutation array used to swap pre with post variables */
  std::unique_ptr<int[]> m_permute;
  /* transition relation */
  BDD m_tr;  
  /* transition relation with m_cube_post abstracted */
  BDD m_tr_nopost;  
  /* BDD cubes with input and post variables */
  BDD m_cube_post;
  BDD m_cube_input;
public:
  /** @brief initialize the enforcabel predecessor
   *  
   * @param manager - the Cudd manager
   * @param transition_relation - the BDD encoding the transition function of the SymbolicModel\n 
   *                              computed with SymbolicModel::compute_gb
   * @param  model - SymbolicModel containing the SymbolicSet for the state and input alphabet 
   **/
  template<class state_type, class input_type>
  EnfPre(const Cudd& manager, 
         const BDD& transition_relation,
         const SymbolicModel<state_type,input_type>& model) : m_tr(transition_relation) {
    /* the permutation array */
    size_t size = manager.ReadSize();
    m_permute = std::unique_ptr<int[]>(new int[size]);
    std::iota(m_permute.get(),m_permute.get()+size,0);
    auto pre_ids = model.get_sym_set_pre().get_bdd_var_ids();
    auto post_ids = model.get_sym_set_post().get_bdd_var_ids();
    for(size_t i=0; i<pre_ids.size(); i++)
      m_permute[pre_ids[i]]=post_ids[i];
    /* create a cube with the input bdd vars */
    m_cube_input = manager.computeCube(model.get_sym_set_input().get_bdd_vars());
    /* create a cube with the post bdd vars */
    m_cube_post = manager.computeCube(model.get_sym_set_post().get_bdd_vars());
    /* copy the transition relation */
    m_tr_nopost=m_tr.ExistAbstract(m_cube_post);
  }
  /** @brief computes the enforcable predecessor of the BDD Z **/
  BDD operator()(BDD Z) const {
    /* project onto state alphabet */
    Z=Z.ExistAbstract(m_cube_post*m_cube_input);
    /* swap variables */
    Z=Z.Permute(m_permute.get());
    /* find the (state, inputs) pairs with a post outside the safe set */
    BDD F = m_tr.AndAbstract(!Z,m_cube_post); 
    /* the remaining (state, input) pairs make up the pre */
    BDD preZ= m_tr_nopost & (!F);
    return preZ;
  }
  const BDD& cube_input() const {
    return m_cube_input;
  }
};


/** 
 *  @brief:  
 *
 **/
inline 
BDD solve_invariance_game(const Cudd& manager, const EnfPre& enf_pre, const BDD& S, bool verbose=true)  {

  BDD Z = manager.bddZero();
  BDD ZZ = manager.bddOne();

  /* as long as not converged */
  size_t i;
  for(i=1; ZZ != Z; i++ ) {
    Z=ZZ;
    ZZ=enf_pre(Z) & S;
    /* print progress */
    if(verbose) {
      std::cout << ".";
      std::flush(std::cout);
      if(!(i%80))
        std::cout << std::endl;
    }
  }
  if(verbose) 
    std::cout << "\nNumber of iterations: " << i << std::endl;
  return Z;
} 


/** 
 *  @brief:  
 *
 **/
inline 
BDD solve_reachability_game(const Cudd& manager, const EnfPre& enf_pre, const BDD& T, bool verbose=true)  {

  BDD Z = manager.bddOne();
  BDD ZZ = manager.bddZero();
  /* the controller */
  BDD C = manager.bddZero();
  /* as long as not converged */
  size_t i;
  for(i=1; ZZ != Z; i++ ) {
    Z=ZZ;
    ZZ=enf_pre(Z) | T;
    /* new (state/input) pairs */
    BDD N = ZZ & (!(C.ExistAbstract(enf_pre.cube_input())));
    /* add new (state/input) pairs to the controller */
    C=C | N;
    /* print progress */
    if(verbose) {
      std::cout << ".";
      std::flush(std::cout);
      if(!(i%80))
        std::cout << std::endl;
    }
  }
  if(verbose) 
    std::cout << "\nNumber of iterations: " << i << std::endl;
  return C;
}

} /* close namespace */
#endif /* FIXEDPOINT_HH_ */
