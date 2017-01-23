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
#include <numeric>

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
  /* a vector of BddIntegerIntervals - one for each dimension */
  std::vector<BddIntegerInterval<abs_type>> m_bdd_interval;
public:
  ///* @cond  EXCLUDE from doxygen */
  ///* default destructor */
  //~SymbolicSet() { };
  ///* default move constructor */
  ////SymbolicSet(SymbolicSet&& other) : m_manager(other.m_manager) {
  ////  *this=std::move(other);
  ////}
  ///* default copy constructor*/
  //SymbolicSet(const SymbolicSet& other)   {
  //  *this=other;
  //}
  ///* default move assignment operator */
  ////SymbolicSet& operator=(SymbolicSet&& other) {
  ////  UniformGrid::operator=(std::move(other));
  ////  m_bdd_interval = std::move(other.m_bdd_interval);
  ////  return *this;
  ////}
  ///* default copy assignment operator */
  //SymbolicSet& operator=(const SymbolicSet& other) {
  //  UniformGrid::operator=(other);
  //  m_bdd_interval = other.m_bdd_interval;
  //  return *this;
  //}
  /* @endcond */
  
  SymbolicSet(const SymbolicSet& other, std::vector<int> dim) :
              UniformGrid(other,dim), m_bdd_interval{} {
    for(const auto& i : dim) 
      m_bdd_interval.push_back(other.m_bdd_interval[i]);
  }
  SymbolicSet(const UniformGrid& grid,
              const std::vector<BddIntegerInterval<abs_type>>& intervals) :
              UniformGrid(grid), m_bdd_interval(intervals) { 
  }
  /* @endcond */

  /** @brief construct SymbolicSet with a Cudd manager **/
  SymbolicSet() : UniformGrid(), m_bdd_interval{} { }
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
              UniformGrid(dim,lb,ub,eta), m_bdd_interval{} {
    for(int i=0; i<m_dim; i++) 
      m_bdd_interval.emplace_back(manager,abs_type{0},m_no_grid_points[i]-1);
  }
  /**
   * @brief provide BDD variable manager and UniformGrid 
   * 
   * @param manager  - BDD variable manager 
   * @param grid     - UnfiormGrid
   **/
  SymbolicSet(const Cudd& manager,
              const UniformGrid& grid) :
              UniformGrid(grid), m_bdd_interval{} {
    for(int i=0; i<m_dim; i++) 
      m_bdd_interval.emplace_back(manager,abs_type{0},m_no_grid_points[i]-1);
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
  SymbolicSet(const SymbolicSet& set1, const SymbolicSet& set2) {
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
    BDD bdd;
    abs_type num;
    for(int k=m_dim-1; k > 0; k--) {
      num=id/m_NN[k];
      id=id%m_NN[k];
      if(k==(m_dim-1))
        bdd = m_bdd_interval[k].int_to_bdd(num);
      else 
        bdd = bdd & m_bdd_interval[k].int_to_bdd(num);
    }
    num=id;
    bdd = bdd & m_bdd_interval[0].int_to_bdd(num);
    return bdd;
  }

  /** @brief get a BDD representation of the grid points  whose IDs are 
   *  an element of the integer hyper-interval [lb; ub] **/
  BDD interval_to_bdd(const Cudd& manager,
                      const std::vector<abs_type>& lb,   
                      const std::vector<abs_type>& ub) const {
    BDD bdd = m_bdd_interval[0].interval_to_bdd(manager,lb[0],ub[0]);
    for(int i=1; i<m_dim; i++) 
      bdd = bdd & m_bdd_interval[i].interval_to_bdd(manager,lb[i],ub[i]);
    return bdd;
  }

  /**
   * @brief obtain a BDD representation of the grid points whose grid point IDs
   * evaluate to true in the lambda expression
   * \verbatim [](const abs_type& i) -> bool \endverbatim
   **/
  template<class F>
  BDD ap_to_bdd(const Cudd& manager, const F& atomic_prop) const {
    BDD result = manager.bddZero();; 
    for(abs_type i=0; i<size(); i++) {
      if(atomic_prop(i)) 
        result = result | id_to_bdd(i);
    }
    return result;
  }


  /** @brief get a vector of grid points that are encoded in the BDD **/
  std::vector<std::vector<double>> bdd_to_grid_points(const Cudd& manager, BDD bdd) {
    if((!get_no_bdd_vars()) || bdd==manager.bddZero())
      return {};
    /* disable reordering (if enabled) */
    Cudd_ReorderingType *method=nullptr;
    if(manager.ReorderingStatus(method))
      manager.AutodynDisable();
    /* get variable ids */
	  auto var_id = get_bdd_var_ids();
	  std::vector<unsigned int> t_ids {};
    for(const auto& v_id : var_id) 
      for(const auto& id : v_id) 
        t_ids.push_back(id);
		/* find the variables in the support of the BDD but outside the SymbolicSet */
		auto support_id = bdd.SupportIndices();
    std::vector<BDD> out{}; 
    for(const auto& id : support_id) {
        if(std::find(t_ids.begin(), t_ids.end(), id)==t_ids.end())
          out.emplace_back(manager.bddVar(id));
    }
    /* remove those variables from the bdd */
    if(out.size()) 
      bdd = bdd.ExistAbstract(manager.computeCube(out));
    /* limit the grid points in the grid */
    for(const auto& interval : m_bdd_interval) 
      bdd = bdd & interval.get_all_elements();
    /* init the vector of grid points to be returned */
    std::vector<double> first(m_dim);
    for(int i=0; i<m_dim; i++)
      first[i]=m_first[i];
    std::vector<std::vector<double>> gp(get_size(bdd),first);
		/* set up iteration to iterate over BDD cubes */
		DdManager* dd = manager.getManager();
	  int *cube;
    CUDD_VALUE_TYPE value;
    DdGen *gen;
    abs_type counter=0;
		/* iterate over BDD cubes */
 		Cudd_ForeachCube(dd,bdd.getNode(),gen,cube,value) {
      abs_type offset=1;
      for(int i=0; i<m_dim; i++) {
        unsigned int no_vars = var_id[i].size();
        for (unsigned int j=0; j<no_vars; j++) {
          if(cube[var_id[i][j]]==1) {
            gp[counter][i]+=(abs_type{1}<<(no_vars-1-j))*m_eta[i];
          }
          /* take care of don't care */
          if(cube[var_id[i][j]]==2) {
            for(abs_type k=0; k<offset; k++) {
              for(int l=0; l<=i; l++) {
                gp[counter+k+offset][l]=gp[counter+k][l];
              }
            }
            for(abs_type k=0; k<offset; k++) {
              gp[counter+k+offset][i]+=(abs_type{1}<<(no_vars-1-j))*m_eta[i];
            }
            offset=(offset<<1);
          }
        }
      }
      counter+=offset;
    }
    /* reactivate reordering if it was enabled */
    if(method!=nullptr)
      manager.AutodynEnable(*method);
    return gp;
  }

  /** @brief projection of the set of grid points onto the dimension in dim+1,
   * e.g., to project onto the 1st  and 2nd dimension, dim = {0,1}  **/
  std::vector<std::vector<double>> projection(const Cudd& manager, const BDD& bdd, std::vector<int> dim) const {
    if(!dim.size())
      return {};
    return SymbolicSet(*this,dim).bdd_to_grid_points(manager,bdd);
  }

  /**
   * @brief restriction of the set of grid points (represented as the BDD bdd) to x
   *
   * Consider the set \f$S\f$ represented by the BDD bdd, being a subset of the
   * n-dimensional hyper-interval 
   * \f$
   *  S\subseteq ( [lb_0,ub_0] \times \cdots [lb_{n-1},ub_{n-1}] ) \cap \eta  \mathbb Z^n
   * \f$.
   * Let \f$x=(x_{j_0},\ldots, x_{j_{m-1}})\f$ with \f$m<n\f$ and \f$p=n-m\f$. In
   * the default setting, the indices \f$j_i\f$ are set to \f$j_i=i\f$. Then the
   * function returns the grid points in the p-dimensional hyper-interval\n\n
   * \f$
   *  S(x):= \{ y\mid (x_{j_0},\ldots, x_{j_m},y_{0},\ldots,y_{p})\in S \}.
   * \f$
   * \n\n
   * Note that \f$S(x) \subseteq [lb_{m},ub_{m}]\times\cdots\times
   * [lb_{n-1},ub_{n-1}] ) \cap \eta  \mathbb Z^n\f$. Optionally, provide the
   * m-dimensional vector domain defining the indices
   * \f$j_0,\ldots,j_{m-1}\f$.
   * For example, if \f$n=4\f$, \f$m=2\f$ and \f$j_0=1,j_1=3\f$, then the
   * function returns the grid points in \n\n
   * \f$
   *  S(x):= \{ (y_0,y_1)\mid (y_0,x_1,y_1,x_3)\in S \}.
   * \f$
   * 
   * @param bdd      the BDD that is used to map 
   * @param x        the element to which the relation is restricted \n
   *                 if domain is not provided, then grid_point_t needs to have grid_point_t::size() method
   * @param domain   Optinally, provide the indices over which the elements of
   *                   x are interpreted
   * @result          vector of grid points  \f$ S(x) \f$
   **/
  template<class grid_point_t>
  std::vector<std::vector<double>> restriction(const Cudd& manager,
                                               const BDD& bdd, 
                                               const grid_point_t& x,  
                                               std::vector<int> domain = {}) const {
    /* compute indices in domain/codomain */
    std::vector<int> codomain {};
    /* fill in default values */
    if(domain.size()==0) {
      domain.resize(x.size());
      std::iota(domain.begin(),domain.end(),0);
      codomain.resize(m_dim-x.size());
      std::iota(codomain.begin(),codomain.end(),x.size());
    } else {
      for(int i=0; i<m_dim; i++) {
        if(std::find(domain.begin(), domain.end(), i) == domain.end())
          codomain.push_back(i);
      }
    }
    /* create SymbolicSet of the domain and codmain */
    SymbolicSet set_dom(*this,domain);
    SymbolicSet set_codom(*this,codomain);
    /* extract bdd with restriction */
    abs_type id = set_dom.xtoi(x);
    BDD restricted = bdd & set_dom.id_to_bdd(id);
    restricted = restricted.ExistAbstract(manager.computeCube(set_dom.get_bdd_vars()));
    /* map restricted BDD to grid points */
    return set_codom.bdd_to_grid_points(manager,restricted);
  }

  /** @brief get number of BDD variables used to represent the SymbolicSet **/
  abs_type get_no_bdd_vars() const {
    abs_type num = 0;
    for(const auto& interval : m_bdd_interval) 
      num+=interval.get_no_bdd_vars();
    return num;
  }

  /** @brief get number of BDD variables used to represent the SymbolicSet **/
  std::vector<BDD> get_bdd_vars() const {
    std::vector<BDD> vars {};
    for(const auto& interval : m_bdd_interval) {
      std::vector<BDD> var = interval.get_bdd_vars();
      for(const auto& p : var) {
        vars.push_back(p);
      }
    }
    return vars;
  }

  /** @brief get number of BDD variables used to represent the SymbolicSet **/
  std::vector<std::vector<unsigned int>> get_bdd_var_ids() const {
    std::vector<std::vector<unsigned int>> var_id {};
    for(const auto& interval : m_bdd_interval) 
      var_id.push_back(interval.get_bdd_var_ids());
    return var_id;
  }


  /** @brief get number grid points represented by the BDD  **/
  abs_type get_size(BDD bdd) const {
    /* limit the grid points in the grid */
    for(const auto& interval : m_bdd_interval) 
      bdd = bdd & interval.get_all_elements();
    return static_cast<abs_type>(bdd.CountMinterm(get_no_bdd_vars()));
  } 

}; /* close class def */
} /* close namespace */
#endif /* SYMBOLICSET_HH_ */
