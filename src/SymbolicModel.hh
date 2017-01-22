/*
 * SymbolicModel.hh
 *
 *  created: Jan 2017
 *   author: Matthias Rungger
 */

/** @file **/
#ifndef SYMBOLICMODEL_HH_
#define SYMBOLICMODEL_HH_

#include <iostream>
#include <vector>


#include "SymbolicSet.hh"

/** @namespace scots **/ 
namespace scots {
/**
 * @class SymbolicModel
 * 
 * @brief Compute the transition function of a symbolic model as BDD using the growth bound. 
 *
 * See 
 * - the manual in <a href="./../../manual/manual.pdf">manual.pdf</a>
 * - http://arxiv.org/abs/1503.03715 for theoretical background 
 *
 **/
template<class state_type, class input_type>
class SymbolicModel {
private:
	/* measurement error bound */
	double* m_z;
  /* print progress to the console (default m_verbose=true) */
  bool m_verbose=true;
	/* SymbolicSet conaining the BDD vars of the pre  */
  const SymbolicSet& m_pre;
	/* SymbolicSet conaining the BDD vars of the inputs  */
  const SymbolicSet& m_input;
	/* SymbolicSet conaining the BDD vars of the post  */
  const SymbolicSet& m_post;

public:
  /* @cond  EXCLUDE from doxygen*/
  /* destructor */
  ~SymbolicModel() {
    delete[] m_z;
  }
  /* deactivate standard constructor */
  SymbolicModel() = delete;
  /* cannot be copied or moved */
  SymbolicModel(SymbolicModel&&) = delete;
  SymbolicModel(const SymbolicModel&) = delete;
  SymbolicModel& operator=(SymbolicModel&&)=delete;
  SymbolicModel& operator=(const SymbolicModel&)=delete;
  /* @endcond */

  /** 
   * @brief Construct SymbolicModel with the SymbolicSet s representing the
   * state alphabet (pre and post) and input alphabet
   * @param  - lambda function with signature  
   *                      \verbatim [] (state_type &x, const input_type &u) ->  void  \endverbatim
   *                      system_post(x,u) provides a numerical approximation of ODE 
   *                      solution at time tau with initial state x and input u \n
   *                      the result is stored in x
   *
   **/
  SymbolicModel(const SymbolicSet& pre,
                const SymbolicSet& input,
                const SymbolicSet& post) :
                m_pre(pre), m_input(input), m_post(post) { 
    m_z = new double[pre.get_dim()] ();
  }

  /** 
   * @brief computes the transition function
   *
   * @param  tf - a BDD encoding the transition function as boolean function over
   *              the BDD var IDs in m_pre, m_input, m_post
   *
   * @param system_post - lambda expression of the form
   *                      \verbatim [] (state_type &x, const input_type &u) ->  void  \endverbatim
   *                      system_post(x,u) provides a numerical approximation of ODE 
   *                      solution at time tau with initial state x and input u \n
   *                      the result is stored in x
   *
   * @param radius_post - lambda expression of the form
   *                      \verbatim [] (state_type &r, const state_type& x, const input_type &u) -> void  \endverbatim
   *                      radius_post(x,u) provides a numerical approximation of
   *                      the growth bound for the cell (with center x, radius  r) and input u\n
   *                      the result is stored in r
   *
   * @result -  number of transitions
   **/
  template<class F1, class F2>
  size_t compute_gb(BDD& tf, F1& system_post, F2& radius_post) {
    return compute_gb(tf, system_post, radius_post, [](const abs_type&) noexcept {return false;});
  }
  /** 
   * @brief computes the transition function
   *
   * @param  tf - a BDD encoding the transition function as boolean function over
   *              the BDD var IDs in m_pre, m_input, m_post
   *
   * @param system_post - lambda expression  as above
   *
   * @param radius_post - lambda expression as above
   *
   * @param avoid    - lambda of the form
   *                      \verbatim [] (abs_type &i) -> bool \endverbatim
   *                      returns true if the abstract state i is in the avoid
   *                      set; otherwise returns false
   * @result -  number of transitions
   **/
  template<class F1, class F2, class F3>
  size_t compute_gb(BDD& tf, F1& system_post, F2& radius_post, F3&& avoid) {
    /* number of cells */
    abs_type N=m_pre.size(); 
    /* number of inputs */
    abs_type M=m_input.size();
    /* state space dimension */
    int dim=m_pre.get_dim();
    /* for display purpose */
    abs_type counter=0;
    /* variables for managing the post */
    std::vector<abs_type> lb(dim);  /* lower-left corner */
    std::vector<abs_type> ub(dim);  /* upper-right corner */
    /* radius of hyper interval containing the attainable set */
    state_type eta;
    state_type r;
    /* state and input variables */
    state_type x;
    input_type u;
    /* for out of bounds check */
    state_type lower_left;
    state_type upper_right;
    /* copy data from m_state_alphabet */
    for(int i=0; i<dim; i++) {
      eta[i]=m_pre.get_eta()[i];
      lower_left[i]=m_pre.get_lower_left()[i];
      upper_right[i]=m_pre.get_upper_right()[i];
    }
    /* the BDD to encode the transition function */
    tf = m_pre.get_zero();
    /* is post of (i,j) out of domain ? */
    bool out_of_domain;
    /* loop over all cells */
    for(abs_type i=0; i<N; i++) {
      BDD bdd_i = m_pre.id_to_bdd(i);
      /* is i an element of the avoid symbols ? */
      if(avoid(i)) {
        continue;
      }
      /* loop over all inputs */
      for(abs_type j=0; j<M; j++) {
        BDD bdd_j = m_input.id_to_bdd(j);
        /* get center x of cell */
        m_pre.itox(i,x);
        /* cell radius (including measurement errors) */
        for(int k=0; k<dim; k++)
          r[k]=eta[k]/2.0+m_z[k];
        /* current input */
        m_input.itox(j,u);
        /* integrate system and radius growth bound */
        /* the result is stored in x and r */
        radius_post(r,x,u);
        system_post(x,u);
        /* determine the cells which intersect with the attainable set: 
         * discrete hyper interval of cell indices 
         * [lb[0]; ub[0]] x .... x [lb[dim-1]; ub[dim-1]]
         * covers attainable set */
        for(int k=0; k<dim; k++) {
          /* check for out of bounds */
          double left = x[k]-r[k]-m_z[k];
          double right = x[k]+r[k]+m_z[k];
          if(left <= lower_left[k]-eta[k]/2.0  || right >= upper_right[k]+eta[k]/2.0) {
            out_of_domain=true;
            break;
          } 
          /* integer coordinate of lower left corner of post */
          lb[k] = static_cast<abs_type>((left-lower_left[k]+eta[k]/2.0)/eta[k]);
          /* integer coordinate of upper right corner of post */
          ub[k] = static_cast<abs_type>((right-lower_left[k]+eta[k]/2.0)/eta[k]);
        }
        if(out_of_domain) {
          out_of_domain=false;
          continue;
        }
        /* compute BDD of post */
        BDD bdd_k = m_post.interval_to_bdd(lb,ub);
        /* add to transition function */
        tf = tf | (bdd_i & bdd_j & bdd_k);
      }
      /* print progress */
      if(m_verbose && ((double)i/(double)N*100)>counter){
        if(counter==0)
          std::cout << "loop: ";
        if((counter%10)==0)
          std::cout << counter;
        else if((counter%2)==0) {
          std::cout << ".";
        }
        counter++;
      }
      std::flush(std::cout); 
    }
    if(m_verbose)
      std::cout << "100" << std::endl;

    /* count number of transitions */
    size_t nvars = m_pre.get_no_bdd_vars() +
                   m_input.get_no_bdd_vars() +
                   m_post.get_no_bdd_vars();
    return tf.CountMinterm(nvars);
  }
//
//  /** @brief TODO
//   *  attainable set associated with cell with center x and input u
//   *
//   *  @param transition_function - the transition function of the abstraction
//   *  @param x - center of cell 
//   *  @param u - input
//   **/
//  void print_post(const state_type& x, const input_type& u) const { }
//
//  /** @brief TODO
//   *  attainable set associated with cell (with center x) and input u
//   *
//   *  @param system_post - lambda function as in compute
//   *  @param radius_post - lambda function as in compute
//   *  @param x - center of cell 
//   *  @param u - input
//   **/
//  template<class F1, class F2>
//  void print_post( ){ }
//
//  /** @brief set the measurement error bound **/
//  void set_measurement_error_bound(const state_type& error_bound) {
//    for(int i=0; i<m_state_alphabet.get_dim(); i++) {
//      m_z[i]=error_bound[i];
//    }
//  }
  /** @brief get measurement error bound **/
  std::vector<double> get_measruement_error_bound() {
    std::vector<double> z;
    for(int i=0; i<m_pre.get_dim(); i++) {
      z.push_back(m_z[i]);
    }
    return z;
  }
  /** @brief activate console output **/
  void verbose_on() {
    m_verbose=true;
  }
  /** @brief deactivate console output **/
  void verbose_off() {
    m_verbose=false;
  }
};

} /* close namespace */
#endif /* SYMBOLICMODEL_HH_ */
