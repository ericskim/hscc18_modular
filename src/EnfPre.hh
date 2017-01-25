/*
 * EnfPre.hh
 *
 *  created: Jan 2017
 *   author: Matthias Rungger
 */

#ifndef ENFPRE_HH_
#define ENFPRE_HH_

#include <iostream>

#include "SymbolicSet.hh"

namespace scots {
/**
 * @class enf_pre
 * 
 * @brief computes of the enforcable predecessor 
 * 
 * Let \f$ F:X\times U \rightrightarrows X\f$ be the transition function and \f$
 * Z\subseteq X\f$, then \n\n
 *
 *  \f$\mathrm{pre}(Z)=\{ (x,u) \in X\times U \mid F(x,u)\neq \emptyset \wedge  F(x,u)\subseteq Z\}\f$
 *
 **/
class enf_pre {
private:
  /* stores the permutation array used to swap pre with post variables */
  std::vector<int> m_permute;

  /* transition relation */
  BDD m_trans_relation;  
  /* transition relation with cubePost_ abstracted */
  BDD m_trans_rel_nopost;  

   /* cubes with input and post variables; used in the existential abstraction  */
  BDD m_cube_post;
  BDD m_cube_input;
  
public:

  /** initialize the enforcabel predecessor object with a <SymbolicModel> containing the
   *  transition relation
   **/
  enf_pre(SymbolicModel *symbolicModel) {
     symbolicModel_=symbolicModel;
    ddmgr_=symbolicModel_->ddmgr_;
     /* the permutation array */
    size_t n=ddmgr_->ReadSize();
    permute_ = new int[n];
    for(size_t i=0; i<n; i++)
      permute_[i]=i;
    for(size_t i=0; i<symbolicModel_->nssVars_; i++)
      permute_[symbolicModel_->preVars_[i]]=symbolicModel_->postVars_[i];
    /* create a cube with the input Vars */
    BDD* vars = new BDD[symbolicModel_->nisVars_];
    for (size_t i=0; i<symbolicModel_->nisVars_; i++)
      vars[i]=ddmgr_->bddVar(symbolicModel_->inpVars_[i]);
    cubeInput_ = ddmgr_->bddComputeCube(vars,NULL,symbolicModel_->nisVars_);
    delete[] vars;
    /* create a cube with the post Vars */
    vars = new BDD[symbolicModel_->nssVars_];
    for (size_t i=0; i<symbolicModel_->nssVars_; i++)
      vars[i]=ddmgr_->bddVar(symbolicModel_->postVars_[i]);   
    cubePost_ = ddmgr_->bddComputeCube(vars,NULL,symbolicModel_->nssVars_);
    delete[] vars;

    /* copy the transition relation */
    R_=symbolicModel_->transitionRelation_;
    RR_=R_.ExistAbstract(cubePost_);
  }
  ~FixedPoint() {
    delete[] permute_;
  }

  BDD operator=(BDD Z)  {
    /* project onto state alphabet */
    Z=Z.ExistAbstract(cube_post*cube_input);
    /* swap variables */
    Z=Z.Permute(permute);
    /* find the (state, inputs) pairs with a post outside the safe set */
    BDD nZ = !Z;
    BDD F = m_trans_relation.AndAbstract(nZ,cube_post); 
    /* the remaining (state, input) pairs make up the pre */
    BDD nF = !F;
    BDD preZ= m_trans_relation_helper.AndAbstract(nF,cube_post);
    return preZ;
  }
}; /* close class def */
} /* close namespace */

#endif /* ENFPRE_HH_ */
