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
  /* @cond  EXCLUDE from doxygen */
  SymbolicSet(const SymbolicSet& other, std::vector<int> dim) :
              UniformGrid(other,dim), m_manager(other.m_manager) {
    for(size_t i=0; i<dim.size(); i++) {
      m_bdd_interval.emplace_back(other.m_bdd_interval[i]);
    }
  }
  /* @endcond */
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

  /** @brief get number of BDD variables used to represent the SymbolicSet **/
  std::vector<std::vector<double>> bdd_to_grid_points(const BDD& bdd) {
    if((!get_no_bdd_vars()) || bdd==get_zero())
      return {};
    /* init return vector of grid points */
    std::vector<std::vector<double>> gp(get_size(bdd),std::vector<double>(m_dim,0));
    /* disable reordering (if enabled) */
    Cudd_ReorderingType *method=nullptr;
    if(m_manager.ReorderingStatus(method))
      m_manager.AutodynDisable();
		
    /* set all BDD variables outside the symbolic set to don't cares */
    std::vector<BDD> vars = get_bdd_vars();
    BDD tmp = vars[0];
    for(auto v : vars) {
      tmp = tmp & (!v);
    }
    /* bdd that is used to extract the grid points */
    BDD bdd_gp = bdd | get_zero();
    /* get variable ids */
    std::vector<std::vector<unsigned int>> var_id(m_dim);
    for(int i=0; i<m_dim; i++) {
      var_id[i]=m_bdd_interval[i].get_bdd_var_ids();
    }
		/* iterate over BDD cubes */
		DdManager* dd = m_manager.getManager();
	  int *cube;
    CUDD_VALUE_TYPE value;
    DdGen *gen;
    abs_type counter=0;
 
 		Cudd_ForeachCube(dd,bdd_gp.getNode(),gen,cube,value) {
      abs_type d=1;
      for(int i=0; i<m_dim; i++) {
        unsigned int no_vars = var_id[i].size();
        for (unsigned int j=0; j<no_vars; j++) {
          if(cube[var_id[i][j]]==1) {
            gp[counter][i]+=(abs_type{1}<<(no_vars-1-j));
          }
          /* take care of don't care */
          if(cube[var_id[i][j]]==2) {
            for(abs_type k=0; k<d; k++) {
              for(int l=0; l<=i; l++) {
                gp[counter+k+d][l]=gp[counter+k][l];
              }
            }
            for(abs_type k=0; k<d; k++) 
                gp[counter+k+d][i]+=(abs_type{1}<<(no_vars-1-j));
            d=(d<<1);/* d=2*d; */
          }
        }
     }
      counter+=d;
    }
    for(auto x : gp) {
      for(auto i : x) 
        std::cout << i << " ";
      std::cout << "\n";
    }

    /* reactivate reordering if it was enabled */
    if(method!=nullptr)
      m_manager.AutodynEnable(*method);
    return gp;
  }
  /**
   * @brief restriction of the set of grid points (represented as the BDD bdd) to x
   * 
   * @param bdd      - the BDD that is used to map 
   * @param x        - the element to which the relation is restricted
   * @param domain   - Optinally, provide the indices over which the elements of
   *                   x are interpreted
   **/
  template<class grid_point_t>
  std::vector<std::vector<double>> restriction(const BDD& bdd, 
                                               const grid_point_t& x,  
                                               std::vector<int> domain = {}) {
    /* init domain/codomain indices */
    std::vector<int> codomain {};
    if(domain.size()==0) {
      for(size_t i=0; i<x.size(); i++) 
        domain.push_back(i);
      for(int i=x.size(); i<m_dim; i++) 
        codomain.push_back(i);
    }
    /* create SymbolicSet of the domain and codmain */
    SymbolicSet set_dom(*this,domain);
    SymbolicSet set_codom(*this,codomain);

    /* extract bdd with restriction */
    abs_type id = set_dom.xtoi(x);
    BDD restriction = bdd & set_dom.id_to_bdd(id);
    restriction = restriction.ExistAbstract(m_manager.computeCube(set_dom.get_bdd_vars()));

    return set_codom.bdd_to_grid_points(restriction);
  }

  /** @brief get number of BDD variables used to represent the SymbolicSet **/
  abs_type get_no_bdd_vars() const {
    abs_type num = 0;
    for(int i=0; i<m_dim; i++) {
      num+=m_bdd_interval[i].get_no_bdd_vars();
    }
    return num;
  }

  /** @brief get number of BDD variables used to represent the SymbolicSet **/
  std::vector<BDD> get_bdd_vars() const {
    std::vector<BDD> vars {};
    for(int i=0; i<m_dim; i++) {
      std::vector<BDD> var = m_bdd_interval[i].get_bdd_vars();
      for(auto p : var) {
        vars.push_back(p);
      }
    }
    return vars;
  }

  /** @brief get number of BDD variables used to represent the SymbolicSet **/
  std::vector<unsigned int> get_bdd_var_ids() const {
    std::vector<unsigned int> var_id {};
    for(int i=0; i<m_dim; i++) {
      std::vector<unsigned int> id = m_bdd_interval[i].get_bdd_var_ids();
      for(auto p : id) {
        var_id.push_back(p);
      }
    }
    return var_id;
  }


  /** @brief get number grid points represented by the BDD  **/
  abs_type get_size(BDD bdd) const {
    return static_cast<abs_type>(bdd.CountMinterm(get_no_bdd_vars()));
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
