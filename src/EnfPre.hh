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

void print_support(const Cudd& mgr, const BDD& x){
  std::vector< unsigned int >  indices = mgr.SupportIndices({x});
  for (size_t i = 0; i < indices.size(); i++){
    std::cout << indices[i] << " ";
  }
  std::cout << std::endl;
}

/** @brief Variables IDs that the BDD x depends on **/
std::vector< unsigned int > get_support(const Cudd& mgr, const BDD& x){
  return mgr.SupportIndices({x});
}

/** @brief Checks whether a BDD formula depends on any of the variables IDs vars. 

If formula is "true" or "false", then automatically dependent. 
**/
bool is_dependent(const Cudd& mgr, const BDD& formula, const BDD& vars){
  if (formula == mgr.bddOne() || formula == mgr.bddZero()){
    return true; 
  }

  std::vector<unsigned int> term_support = get_support(mgr, formula);
  std::sort(term_support.begin(), term_support.end());
  std::vector<unsigned int> var_support = get_support(mgr, vars);
  std::sort(var_support.begin(), var_support.end());
  std::vector<unsigned int>::iterator iter;
  std::vector<unsigned int> intersection_IDs(var_support.size() + term_support.size());
  iter = std::set_intersection(term_support.begin(), term_support.end(), 
                               var_support.begin(), var_support.end(), 
                               intersection_IDs.begin());
  // Check if the intersection of the support is nonempty
  if (iter - intersection_IDs.begin() > 0){
    return true;
  }
  return false; 
}

/**
@brief Applies an existential quantification over a conjunction of BDD formulas. Takes into account variable dependencies. 

@param elim_vars [in] -  variables to eliminate
@param remaining_formulas [in] - copy of  initial formulas in the conjunction
**/
BDD exists_over_conjunction(const Cudd& mgr, std::list<BDD> elim_vars, std::vector<BDD> remaining_formulas) {
  /* TODO  */
  BDD prev = mgr.bddOne();
  std::cout << "No cost function is currently used for variable ordering" << std::endl;
  std::cout << "Formulas in conjunction: " << remaining_formulas.size() << std::endl;
  while(elim_vars.size() > 0){
    /** Identify which variable to eliminate **/
    std::cout << "Remaining Variables: " << elim_vars.size() << std::endl;
    print_support(mgr, prev);
    BDD to_elim = elim_vars.front();
    elim_vars.pop_front();
    // std::cout << "Abstracting out: " << std::endl;
    // print_support(mgr, to_elim);

    /** Identify which formulas are independent **/

    auto is_dependent_on_to_elim = std::bind(is_dependent, mgr, std::placeholders::_1, to_elim); // create predicate
    std::vector<BDD>::iterator bound = std::partition(remaining_formulas.begin(), remaining_formulas.end(), is_dependent_on_to_elim); 

    /** Compute conjunction and quantify variables out**/
    for(auto i=remaining_formulas.begin(); i != bound; i++){
      prev &= *i;
    }
    prev = prev.ExistAbstract(to_elim);

    // Get rid of formulas used in the conjunction 
    bound = std::remove_if(remaining_formulas.begin(), remaining_formulas.end(), is_dependent_on_to_elim);
    remaining_formulas.resize(bound - remaining_formulas.begin());

  }
  std::cout << "No remaining variables" << std::endl;
  print_support(mgr, prev);
  return prev;
}

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
  /** @brief initialize the enforcable predecessor
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

/**
@brief Predecessor operation that attempts to minimize the size of the intermediate BDD.

Unlike EnfPre, it takes all of the sets and the system transition relation as decomposed sets. 
**/
class DecomposedPredecessor{
private:

 

protected: 
  /* stores the permutation array used to swap pre with post variables */
  std::unique_ptr<int[]> m_permute;

  const std::vector<BDD> relation_conjunction;
  std::vector<BDD> m_cubes_input;
  std::vector<BDD> m_cubes_latent;
  std::vector<BDD> m_cubes_post;
  std::list<BDD> pre_elim_vars; 
  Cudd mgr;
public: 
  DecomposedPredecessor(const Cudd& manager,
                        const std::vector<BDD>& relations,
                        const std::vector<SymbolicSet>& pre_sets,
                        const std::vector<SymbolicSet>& control_sets,
                        const std::vector<SymbolicSet>& post_sets,
                        const std::vector<SymbolicSet>& latent_sets): relation_conjunction(relations){
    mgr = manager; 
    std::cout << "Permutation Array" << std::endl;
    /* Permutation array used to swap pre and post states */
    size_t size = mgr.ReadSize(); 
    m_permute = std::unique_ptr<int[]>(new int[size]);
    std::iota(m_permute.get(),m_permute.get()+size,0);
    for (size_t i = 0; i < pre_sets.size(); i++){
      auto pre_ids = pre_sets[i].get_bdd_var_ids();
      auto post_ids = post_sets[i].get_bdd_var_ids();
      for(size_t i=0; i<pre_ids.size(); i++)
        m_permute[pre_ids[i]]=post_ids[i];
    }

    m_cubes_post.resize(post_sets.size());
    m_cubes_input.resize(control_sets.size());
    m_cubes_latent.resize(latent_sets.size());

    /* Compute BDD cubes of different variable domains */
    std::cout << "BDD state Cubes" << std::endl;
    for (size_t i = 0; i < pre_sets.size(); i++){
      std::cout << i << std::endl;
      m_cubes_post[i] = post_sets[i].get_cube(mgr);
      pre_elim_vars.push_back(m_cubes_post[i]); 
    }
    std::cout << "BDD control Cubes" << std::endl;
    for (size_t i = 0; i < control_sets.size(); i++){
      m_cubes_input[i] = control_sets[i].get_cube(mgr);
    }
    std::cout << "BDD latent Cubes" << std::endl;
    for (size_t i = 0; i < latent_sets.size(); i++){
      m_cubes_latent[i] = latent_sets[i].get_cube(mgr);
      pre_elim_vars.push_back(m_cubes_latent[i]);
    }

  }

  /** @brief Computes a robust enforcing pre_state-input pairs that
      ensure the post_state is within Z
  **/
  BDD operator() (BDD Z ) const {
    /* Change Z from a formula over predecessor states to one over post states */
    Z=Z.Permute(m_permute.get());
    std::vector<BDD> conj_formulas(relation_conjunction);
    conj_formulas.push_back(!Z);
    return !exists_over_conjunction(mgr, pre_elim_vars, conj_formulas);
  }

  /** @brief Compute non-blocking state-input pairs **/
  BDD nonblocking() const{
    return exists_over_conjunction(mgr, pre_elim_vars, relation_conjunction);
  } 
};

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
