/*
 *  EnfPre.hh
 *
 *  created: Jan 2017
 *   author: Matthias Rungger
 */

/** @file **/

#ifndef ENFPRE_HH_
#define ENFPRE_HH_

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
protected:
  /* stores the permutation array used to swap pre with post variables */
  std::unique_ptr<int[]> m_permute;
  /* transition relation */
  const BDD m_tr;
  /* transition relation with m_cube_post abstracted */
  BDD m_tr_nopost;  
  /* BDD cubes with input and post variables */
  BDD m_cube_post;
  BDD m_cube_input;
  Cudd mgr;
public:
  /** @brief initialize the enforcabel predecessor
   *  
   * @param manager - the Cudd manager
   * @param transition_relation - the BDD encoding the transition function of the SymbolicModel\n 
   *                              computed with SymbolicModel::compute_gb
   * @param  model - SymbolicModel containing the SymbolicSet for the state and input alphabet 
   **/
  EnfPre(const Cudd& manager, 
         const BDD& transition_relation,
         const SymbolicSet& pre_set,
         const SymbolicSet& control_set,
         const SymbolicSet& post_set) : m_tr(transition_relation){
    mgr = manager; 
    /* the permutation array */
    size_t size = manager.ReadSize();
    m_permute = std::unique_ptr<int[]>(new int[size]);
    std::iota(m_permute.get(),m_permute.get()+size,0);
    auto pre_ids = pre_set.get_bdd_var_ids();
    auto post_ids = post_set.get_bdd_var_ids();
    for(size_t i=0; i<pre_ids.size(); i++)
      m_permute[pre_ids[i]]=post_ids[i];
    /* create a cube with the input bdd vars */
    m_cube_input = control_set.get_cube(manager);
    /* create a cube with the post bdd vars */
    m_cube_post = post_set.get_cube(manager);
    /* non blocking states-input pairs */
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
  
  /** @brief: small function to output progess of an iteration to the terminal **/
  inline void print_progress(int i) {
    std::cout << ".";
    std::flush(std::cout);
    if(!(i%40)) {
      std::cout << "\r";
      std::cout << "                                        ";
      std::cout << "\r";
    }
  }

}; // close EnfPre 



class InterconnectedEnfPre: public EnfPre {
private:
  BDD inter;
  std::vector<BDD> systems;
  BDD m_cube_exog; 
public:
  /** @brief initialize the enforcabel predecessor
   *  
   * @param manager - the Cudd manager
   * @param transition_relation - the BDD encoding the transition function of the SymbolicModel\n 
   *                              computed with SymbolicModel::compute_gb
   * @param  model - SymbolicModel containing the SymbolicSet for the state and input alphabet 
   **/
  InterconnectedEnfPre(const Cudd& manager, 
         const BDD& transition_relation,
         const SymbolicSet& pre_set,
         const SymbolicSet& control_set,
         const SymbolicSet& post_set,
         const SymbolicSet& exog_set,
         const BDD interconnection,
         const std::vector<BDD>& sys_vec) : EnfPre(manager, transition_relation, pre_set, control_set, post_set), 
                                       inter(interconnection),
                                       systems(sys_vec){
    m_cube_exog = exog_set.get_cube(manager);
    m_tr_nopost = m_tr_nopost.ExistAbstract(m_cube_exog); 
  }

  /** @brief computes the enforcable predecessor of the BDD Z **/
  BDD operator()(BDD Z) const {
    // TODO there might be a bug where it's permissible to go outside the interconnection's domain
    /* project onto pre state alphabet */
    Z=Z.ExistAbstract(m_cube_post*m_cube_input*m_cube_exog);
    /* swap pre variables to post. Z is now in post domain */
    Z=Z.Permute(m_permute.get());
    /* find the set of (state, exog, inputs) tuples F with a post intersecting the unsafe set */
    BDD F = m_tr.AndAbstract(!Z,m_cube_post);
    /* get rid of transitions that are inconsistent with the interconnection relation */
    F &= inter;
    F = F.ExistAbstract(m_cube_exog*m_cube_post);
    /* the remaining (state, input) pairs make up the pre */
    return  m_tr_nopost & (!F);
  }
}; // close InterconnectedEnfPre 

//inline 
//BDD solve_invariance_game(const Cudd& manager, const EnfPre& enf_pre, const BDD& S, bool verbose=true)  {
//
//  BDD Z = manager.bddZero();
//  BDD ZZ = manager.bddOne();
//
//  /* as long as not converged */
//  size_t i;
//  for(i=1; ZZ != Z; i++ ) {
//    Z=ZZ;
//    ZZ=enf_pre(Z) & S;
//    /* print progress */
//    if(verbose) {
//      std::cout << ".";
//      std::flush(std::cout);
//      if(!(i%80))
//        std::cout << std::endl;
//    }
//  }
//  if(verbose) 
//    std::cout << "\nNumber of iterations: " << i << std::endl;
//  return Z;
//} 

} /* close namespace */
#endif /* ENFPRE_HH_ */
