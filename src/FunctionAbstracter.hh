/*
 * SymbolicModel.hh
 *
 *  created: Sep 2017
 *   author: Eric S. Kim
 *           
 */

/** @file **/
#ifndef FunctionAbstracter_HH_
#define FunctionAbstracter_HH_

#include <iostream>
#include <vector>
#include <functional>
#include <string>
#include <map>

#include "SymbolicSet.hh"
//#include "Dependency.hh"

/** @namespace scots **/ 
namespace scots {

/**
@brief Keeps track of the dependencies between input/output variables of a given function.
**/
class FunctionDependency{
protected:
  std::vector<SymbolicSet> i_spaces;
  std::vector<SymbolicSet> o_spaces;
  std::vector<int> i_dims, o_dims;

  SymbolicSet i_product, o_product;
  /*@brief Product space dimensions */
  int i_dim, o_dim;
  /*@brief Vector of sets of state variables located on the right hand side of update equation*/
  std::vector<std::vector<int> > depends;

private: 

  /*@brief Maps a output space set and coordinate to an index in pre_deps*/
  //std::map<IntegerInterval<abs_type>, size_t> interval_to_index;
  std::vector<std::pair<IntegerInterval<abs_type>, int> > i_interval_to_index;
  std::vector<std::pair<IntegerInterval<abs_type>, int> > o_interval_to_index;

  int get_index(std::vector<std::pair<IntegerInterval<abs_type>, int> > interval_to_index, 
                IntegerInterval<abs_type> set){
    for (size_t i = 0; i < interval_to_index.size(); i++){
      if (interval_to_index[i].first == set)
        return interval_to_index[i].second;
    }
    /*Index not found */
    return -1;
  }

public:
  /**
  @brief FunctionDependency constructor.
  @param [in] isets - The function's input set is a cartesian product of symbolic sets in isets
  @param [in] osets - The function's output set is a cartesian product of symbolic sets in osets
  **/
  FunctionDependency(){};
  FunctionDependency(std::vector<SymbolicSet> isets, std::vector<SymbolicSet> osets): 
                  i_spaces(isets), o_spaces(osets){
    i_dims.resize(i_spaces.size());
    o_dims.resize(o_spaces.size());
    for (size_t i = 0; i < i_spaces.size(); i++){
      i_dims[i] = i_spaces[i].get_dim();
    }
    for (size_t i = 0; i < o_spaces.size(); i++){
      o_dims[i] = o_spaces[i].get_dim();
    }

    i_product = i_spaces[0];
    i_dim = i_spaces[0].get_dim();
    for (size_t i = 1; i < i_spaces.size(); i++){
      i_product = SymbolicSet(i_product, i_spaces[i]);
      i_dim += i_spaces[i].get_dim(); 
    }
    if (i_product.get_dim() != i_dim)
      throw std::runtime_error("\nscots::FunctionDependency Input dimensions don't match");
    

    o_product = o_spaces[0];
    o_dim = o_spaces[0].get_dim();
    for (size_t i = 1; i < o_spaces.size(); i++){
      o_product = SymbolicSet(o_product, o_spaces[i]);
      o_dim += o_spaces[i].get_dim();
    }
    if (o_product.get_dim() != o_dim)
      throw std::runtime_error("\nscots::FunctionDependency Output dimensions don't match");
    

    /* Associate each IntegerInterval in input product set with a unique index */
    for (int i = 0; i < i_dim; i++){
      i_interval_to_index.push_back(std::make_pair(i_product[i], i));
    }
    for (int i = 0; i < o_dim; i++){
      o_interval_to_index.push_back(std::make_pair(o_product[i], i));
    }

    depends.resize(o_dim);
    for (int i = 0; i < o_dim; i++){
      depends[i].clear();
    }

  }

  ~FunctionDependency(){};

 /** @brief Set function dependency
  If a dependency has already been set for the specific output then it will be cleared if written.
 **/
  void set_dependency(IntegerInterval<abs_type>& oslice, 
                      std::vector<IntegerInterval<abs_type> > pre_deps){
    // TODO Error handling 
    int o_index = get_index(o_interval_to_index, oslice);
    depends[o_index].clear();

    int p_index;
    scots::SymbolicSet iset; 
    for (size_t i = 0; i < pre_deps.size(); i++){
      p_index = get_index(i_interval_to_index, pre_deps[i]);
      depends[o_index].push_back(p_index);
    }

  }

  std::vector<std::vector<int> > get_dependency() const{
    return depends;
  };

  /** @brief Get vector of input dimensions for which output dimension 
             i depends
  **/
  std::vector<int> get_dependency(int i) const{
    return depends[i];
  }

  SymbolicSet get_output_product() const {
    return o_product;
  }

  SymbolicSet get_input_product() const {
    return i_product;
  }

  /** @brief Print out dependencies for the product set **/
  friend std::ostream &operator<< (std::ostream &os, const FunctionDependency & dep){
    for(int i = 0; i < dep.o_dim; i++){
      os << "\nOutput index " << i << " Dependencies: ";
      for (auto &j: dep.depends[i]){
        os << "  " << j;
      }
    }
    os << "\n";

    return os;
  } 

};

/**
@class FunctionAbstracter
@brief 
**/
template<class concreteInput, class concreteOutput>
class FunctionAbstracter{
private:
  /* SymbolicSet containing the BDD vars of the output space */
  SymbolicSet m_outSpace;
  /* SymbolicSet containing the BDD vars of the input space */
  SymbolicSet m_inSpace;
  /* Explicit dependency graph*/
  FunctionDependency dep;
  /* Overapproximates concrete function */
  std::function<void(const concreteInput&, const concreteInput&, concreteOutput &, concreteOutput &)> overApprox;
  
  /**
  *  @brief Lifts the i-th gridpoint in lower dimension space "small" to a full dimensional gridpoint
  **/
  template <class T>
  inline T lifted_input(abs_type i, 
                        const SymbolicSet & full, 
                        const SymbolicSet & small, 
                        std::vector<int> dep,
                        std::string corner = "center"){
    T x;
    static std::vector<double> proj_x;

    /* Dummy values for full dimensional input that don't affect output */
    for (int k=0; k < full.get_dim(); k++)
      x[k] = (full.get_center())[k];

    proj_x.resize(small.get_dim());
    if (corner == "ur"){
      small.i_to_ur(i,proj_x);
    } else if (corner == "ll"){
      small.i_to_ll(i,proj_x);
    } else{
      small.itox(i,proj_x);
    }
    /* Substitute in correct values for important inputs */
    for (int k=0; k < small.get_dim(); k++)
      x[dep[k]] = proj_x[k];
    return x;
  }

  /**
  *  @brief Helper function to find integer coordinates along the post_dim-th coordinate
  * @param[out] post_lb
  * @param[out] post_ub
  * @return  false if out of bounds, true otherwise
  **/
  bool post_interval_bounds(int odim, 
                            concreteOutput ll, concreteOutput ur, 
                            concreteOutput ll_bound, concreteOutput ur_bound, 
                            abs_type &post_lb, abs_type &post_ub){
    std::vector<double> eta = m_outSpace.get_eta();
    if(ll[odim] <= (ll_bound[odim]-eta[odim]/2.0)  || ur[odim] >= (ur_bound[odim]+eta[odim]/2.0)){
      std::cout << "Post dim " << odim << std::endl;
      std::cout << "LL " << ll[odim] << " " << (ll_bound[odim]-eta[odim]/2.0) << std::endl;
      std::cout << "UR " << ur[odim] << " " << (ur_bound[odim]+eta[odim]/2.0) << std::endl;
      return false;
    }
    /* integer coordinate of lower left corner of output dimension odim */
    post_lb = static_cast<abs_type>((ll[odim]-ll_bound[odim]+eta[odim]/2.0)/eta[odim]);
    /* integer coordinate of upper right corner of output dimension odim */
    post_ub = static_cast<abs_type>((ur[odim]-ll_bound[odim]+eta[odim]/2.0)/eta[odim]);
    if (post_ub < post_lb){
      throw std::runtime_error("scots::FunctionAbstracter Unexpected interval index");
    }
    return true;
  }

public:
  ~FunctionAbstracter() {};
  // /* deactivate standard constructor */
  FunctionAbstracter() {};
  // /* cannot be copied or moved */
  // FunctionAbstracter(FunctionAbstracter&&) = delete;
  // FunctionAbstracter(const FunctionAbstracter&) = delete;
  // FunctionAbstracter& operator=(FunctionAbstracter&&)=delete;
  // FunctionAbstracter& operator=(const FunctionAbstracter&)=delete;

  FunctionAbstracter(const FunctionDependency d,
                      const std::function<void(const concreteInput&, const concreteInput&, concreteOutput &, concreteOutput &)> oa): dep(d) {
      m_outSpace = dep.get_output_product();
      m_inSpace = dep.get_input_product();
      overApprox = oa;
  }

  BDD compute_abstraction(const Cudd& mgr){
    const int odims = m_outSpace.get_dim();
    const int idims = m_inSpace.get_dim();
    
    /* variables for managing the low dimensional output */
    abs_type post_lb, post_ub;

    concreteInput input_ll, input_ur;

    /* for out of bounds check on output */
    concreteOutput set_ll, set_ur;
    concreteOutput overapprox_ll, overapprox_ur;
    /* copy data from m_state_alphabet */
    for(int i=0; i<odims; i++) {
      //eta[i]=m_inSpace.get_eta()[i];
      set_ll[i]=m_outSpace.get_lower_left()[i];
      set_ur[i]=m_outSpace.get_upper_right()[i];
    }

    /* the BDD to encode the approximation*/
    BDD approx = mgr.bddOne();
    
    if (odims <= 0){
      throw std::runtime_error("scots::FunctionAbstracter Output space should not be zero dimensional");
    }

    for(int post_dim=0; post_dim<odims; post_dim++){
      BDD coord_approx = mgr.bddZero();
      /*Iterate over input space and compute overapproximations to output*/
      SymbolicSet dep_set = SymbolicSet(m_inSpace, dep.get_dependency(post_dim));
      SymbolicSet post_dim_slice  = SymbolicSet(m_outSpace, {post_dim});

      abs_type N = dep_set.size();
      for(abs_type i=0; i<N; i++) {
        BDD bdd_i = dep_set.id_to_bdd(i);
        input_ll = lifted_input<concreteInput>(i, m_inSpace, dep_set, dep.get_dependency(post_dim), "ll");
        input_ur = lifted_input<concreteInput>(i, m_inSpace, dep_set, dep.get_dependency(post_dim), "ur");

        /* Compute concrete values for post_lower and post_upper */
        overApprox(input_ll, input_ur, overapprox_ll, overapprox_ur);

        bool verbose = true;
        /* Check for out of bounds errors along post_dim*/
        if (!post_interval_bounds(post_dim, overapprox_ll, overapprox_ur, set_ll, set_ur, post_lb, post_ub)){
          std::cout << "Postdim: "<< post_dim << "  Dependency Index: " << i << std::endl;
          std::cout << "In LL: " << input_ll[0] << " " << input_ll[1] << " " << input_ll[2] << " " << input_ll[3] << std::endl; 
          std::cout << "In UR: " << input_ur[0] << " " << input_ur[1] << " " << input_ur[2] << " " << input_ur[3] << std::endl;
          std::cout << "Out LL: " << overapprox_ll[0] << " " << overapprox_ll[1] << std::endl; 
          std::cout << "Out UR: " << overapprox_ur[0] << " " << overapprox_ur[1] << std::endl;
          std::cout << "OUT OF REGION\n";
          std::cout << "Symbolic Indices: " << post_dim << " " << post_lb << " " << post_ub << std::endl << std::endl; 
          continue;
        }

        /* Compute BDD of the post_dim component of the function output */
        BDD bdd_post_component = post_dim_slice.interval_to_bdd(mgr,post_lb,post_ub);


        /* Add transition to current post coordinate */
        coord_approx = coord_approx | (bdd_i & bdd_post_component);
      }
      
      /* Impose coordinate constraint on post_dim*/
      approx &= coord_approx;

    } // end for along output coordinates

    return approx;
  }


}; // class FunctionAbstracter  



} // namespace scots 

#endif /* FunctionAbstracter_HH_ */