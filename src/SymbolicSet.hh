/*
 * SymbolicSet.hh
 *
 *  created: Jan 2016
 *   author: Matthias Rungger
 */

/** @file **/

#ifndef SYMBOLICSET_HH_
#define SYMBOLICSET_HH_


#include <vector>
#include <iostream>
#include <sstream>

#include "UniformGrid.hh"
#include "BddIntegerInterval.hh"

/* cudd library */
#include "cuddObj.hh"

namespace scots {

/** 
 *  @class SymbolicSet 
 *
 *  @brief The UniformGrid class with the capability to represent sets of grid
 *  points as BDDs. See UniformGrid for details.
 *
 **/

class SymbolicSet : public UniformGrid {
private:
  /* reference to BDD manager */
  const Cudd& m_manager;
  /* a vector of BddIntegerIntervals - one for each dimension */
  std::vector<BddIntegerInterval<abs_type>> m_bdd_interval;
public:
  /**
   * @brief provide BDD variable manager and  UniformGrid parameters 
   * 
   * @param manager  - BDD variable manager 
   * @param dim      - dimension of the real space
   * @param lb       - lower-left corner of the hyper-interval confining the uniform grid
   * @param ub       - upper-right corner of the hyper-interval confining the uniform grid
   * @param eta      - grid point distances
   **/
  template<class grid_point_t>
  SymbolicSet(const Cudd& manager,
              const int dim,
              const grid_point_t& lb,
              const grid_point_t& ub,
              const grid_point_t& eta) :
              UniformGrid(dim,lb,ub,eta), m_manager(manager) {
    for(int i=0; i<m_dim; i++) {
      m_bdd_interval.emplace_back(BddIntegerInterval<abs_type>(m_manager,abs_type(0),m_no_grid_points[i]-1));
    }
  }
  /**
   * @brief provide BDD variable manager and UniformGrid 
   * 
   * @param manager  - BDD variable manager 
   * @param grid     - UnfiormGrid
   **/
  SymbolicSet(const Cudd& manager,
              const UniformGrid& grid) :
              UniformGrid(grid), m_manager(manager) {
    for(int i=0; i<m_dim; i++) {
      m_bdd_interval.emplace_back(BddIntegerInterval<abs_type>(m_manager,abs_type(0),m_no_grid_points[i]-1));
    }
  }
  /**
   * @brief construct product of two SymbolicSet s
   * 
   * The instantiated SymbolicSet represents the Cartesian product of the
   * SymbolicSet set1 and the SymbolicSet set2
   * 
   * @param set1  - SymbolicSet
   * @param set2  - SymbolicSet
   **/
  SymbolicSet(const SymbolicSet& set1,
              const SymbolicSet& set2) : m_manager(set1.m_manager) {
    m_dim = set1.m_dim + set2.m_dim;
    m_eta = new double[m_dim];
    m_first = new double[m_dim];
    m_no_grid_points = new abs_type[m_dim];
    m_NN = new abs_type[m_dim];
    for(int i=0; i<set1.m_dim; i++) {
      m_eta[i] = set1.m_eta[i];
      m_first[i]  = set1.m_first[i];
      m_no_grid_points[i]  = set1.m_no_grid_points[i];
      m_bdd_interval.push_back(set1.m_bdd_interval[i]);
    }
    for(int j=set1.m_dim, i=0; i<set2.m_dim; i++) {
      m_eta[j+i] = set2.m_eta[i];
      m_first[j+i]  = set2.m_first[i];
      m_no_grid_points[j+i]  = set2.m_no_grid_points[i];
      m_bdd_interval.push_back(set2.m_bdd_interval[i]);
    }
    calc_nn();
  }

  /** @brief function to obtain a BDD representation of the grid point id **/
  BDD id_to_bdd(abs_type id) const {
    BDD bdd = m_manager.bddOne();
    abs_type num;
    for(int k=m_dim-1; k > 0; k--) {
      num=id/m_NN[k];
      id=id%m_NN[k];
      bdd = bdd & m_bdd_interval[k].int_to_bdd(num);
    }
    num=id;
    bdd = bdd & m_bdd_interval[0].int_to_bdd(num);
    return bdd;
  }

  /** @brief function to obtain a BDD representation of the grid points
   *  whose IDs are an element of the integer hyper-interval [lb; ub] **/
  BDD interval_to_bdd(const std::vector<abs_type>& lb,   
                      const std::vector<abs_type>& ub) const {
    BDD bdd = m_manager.bddOne();
    for(int i=0; i<m_dim; i++) {
      bdd = bdd & m_bdd_interval[i].interval_to_bdd(lb[i],ub[i]);
    }
    return bdd;
  }

  /** @brief print information about the symbolic set **/
  void print_info(int verbose=0) const {
    UniformGrid::print_info();
    std::cout << "Number of BDD variables per dimension ";
    for(int i=0; i<m_dim; i++) {
      std::cout << m_bdd_interval[i].get_no_bdd_vars() << " ";
    }
    if(verbose) {
      std::cout << "\n";
      for(int i=0; i<m_dim; i++) {
        std::cout << "Dim " << i+1 << ": ";
        m_bdd_interval[i].print_bdd_IDs();
      }
    }
    std::cout << "\n";
  }

  /** @brief get number of BDD variables used to represent the SymbolicSet **/
  abs_type get_no_bdd_vars() const {
    abs_type num = 0;
    for(int i=0; i<m_dim; i++) {
      num+=m_bdd_interval[i].get_no_bdd_vars();
    }
    return num;
  }
  BDD get_zero() const {
    return m_manager.bddZero();
  }
  BDD get_one() const {
    return m_manager.bddOne();
  }

}; /* close class def */
} /* close namespace */
#endif /* SYMBOLICSET_HH_ */
