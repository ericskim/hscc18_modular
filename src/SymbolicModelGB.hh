/*
 * SymbolicModelGB.hh
 *
 *  created: Jan 2017
 *   author: Matthias Rungger
 */

/** @file **/
#ifndef SYMBOLICMODELGB_HH_
#define SYMBOLICMODELGB_HH_

#include <iostream>
#include <cstring>
#include <memory>
#include <vector>


#include "UniformGrid.hh"
#include "IndexSet.hh"

/** @namespace scots **/ 
namespace scots {
/**
 * @class SymbolicModelGB
 * 
 * @brief TODO
 *
 * See 
 * - the manual in <a href="./../../manual/manual.pdf">manual.pdf</a>
 * - http://arxiv.org/abs/1503.03715 for theoretical background 
 *
 **/
template<class state_type, class input_type>
class SymbolicModelGB {
private:
	/* measurement error bound */
	double* m_z;
	/* UniformGrid for the state alphabet  */
  const UniformGrid& m_state_alphabet;
	/* UniformGrid for the input alphabet  */
  const UniformGrid& m_input_alphabet;
  /* print progress to the console (default m_verbose=true) */
  bool m_verbose=true;
	/* IndexSet conaining the BDD vars for the pre  */
  const IndexSet& m_pre;
	/* IndexSet conaining the BDD vars for the inputs */
  const IndexSet& m_input;
	/* IndexSet conaining the BDD vars for the post  */
  const IndexSet& m_post;

public:
  /* @cond  EXCLUDE from doxygen*/
  /* destructor */
  ~SymbolicModelGB() {
    delete[] m_z;
  }
  /* deactivate standard constructor */
  SymbolicModelGB() = delete;
  /* cannot be copied or moved */
  SymbolicModelGB(SymbolicModelGB&&) = delete;
  SymbolicModelGB(const SymbolicModelGB&) = delete;
  SymbolicModelGB& operator=(SymbolicModelGB&&)=delete;
  SymbolicModelGB& operator=(const SymbolicModelGB&)=delete;
  /* @endcond */

  /** @brief TODO **/
  SymbolicModelGB(const UniformGrid& state_alphabet,
                const UniformGrid& input_alphabet,
                const IndexSet& bdd_pre,
                const IndexSet& bdd_input,
                const IndexSet& bdd_post) :
                m_state_alphabet(state_alphabet), 
                m_input_alphabet(input_alphabet),
                m_pre(bdd_pre),
                m_input(bdd_input),
                m_post(bdd_post) {

    /* TODO check manager  */
    if((m_pre.get_manager()!=m_input.get_manager()) ||
       (m_input.get_manager()!=m_post.get_manager())) {
      /* error std exception ... */
    }

    m_z = new double[state_alphabet.get_dim()];
    for(int i=0; i<state_alphabet.get_dim(); i++) {
      m_z[i]=0;
    }
  }

  /** 
   * @brief TODO
   *
   * @param system_post - lambda function with signature  
   *                      \verbatim [] (state_type &x, const input_type &u) ->  void  \endverbatim
   *                      system_post(x,u) provides a numerical approximation of ODE 
   *                      solution at time tau with initial state x and input u \n
   *                      the result is stored in x
   *
   * @param radius_post - lambda function with signature
   *                      \verbatim [] (state_type &r, const state_type& x, const input_type &u) -> void  \endverbatim
   *                      radius_post(x,u) provides a numerical approximation of
   *                      the growth bound for the cell (with center x, radius  r) and input u\n
   *                      the result is stored in r
   *
   * 
   **/
  template<class F1, class F2>
  BDD compute(F1& system_post, F2& radius_post) {
    return compute(system_post, radius_post, [](const abs_type&) noexcept {return false;});
  }
  template<class F1, class F2, class F3>
  BDD compute(F1& system_post, F2& radius_post, F3&& overflow) {
    /* number of cells */
    abs_type N=m_state_alphabet.size(); 
    /* number of inputs */
    abs_type M=m_input_alphabet.size();
    /* state space dimension */
    int dim=m_state_alphabet.get_dim();
    /* for display purpose */
    abs_type counter=0;
    /* some grid information */
    std::vector<abs_type> NN=m_state_alphabet.get_nn();
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
      eta[i]=m_state_alphabet.get_eta()[i];
      lower_left[i]=m_state_alphabet.get_lower_left()[i];
      upper_right[i]=m_state_alphabet.get_upper_right()[i];
    }
    
    BDD tf = m_pre.get_zero();
    /* is post of (i,j) out of domain ? */
    bool out_of_domain;
    /* loop over all cells */
    for(abs_type i=0; i<N; i++) {
      BDD i_bdd = m_pre.id_to_bdd(i);
      /* loop over all inputs */
      for(abs_type j=0; j<M; j++) {
        /* is i an element of the overflow symbols ? */
        if(!j & overflow(i)) {
          break;
        }
        BDD j_bdd = m_input.id_to_bdd(j);
        /* get center x of cell */
        m_state_alphabet.itox(i,x);
        /* cell radius (including measurement errors) */
        for(int k=0; k<dim; k++)
          r[k]=eta[k]/2.0+m_z[k];
        /* current input */
        m_input_alphabet.itox(j,u);
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
        BDD k_bdd = m_post.interval_to_bdd(lb,ub);
        /* add to transition function */
        tf = tf | (i_bdd & j_bdd & k_bdd);
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

    return tf;
  }

  /** @brief TODO
   *  attainable set associated with cell with center x and input u
   *
   *  @param transition_function - the transition function of the abstraction
   *  @param x - center of cell 
   *  @param u - input
   **/
  void print_post(const state_type& x, const input_type& u) const { }

  /** @brief TODO
   *  attainable set associated with cell (with center x) and input u
   *
   *  @param system_post - lambda function as in compute
   *  @param radius_post - lambda function as in compute
   *  @param x - center of cell 
   *  @param u - input
   **/
  template<class F1, class F2>
  void print_post( ){ }

  /** @brief set the measurement error bound **/
  void set_measurement_error_bound(const state_type& error_bound) {
    for(int i=0; i<m_state_alphabet.get_dim(); i++) {
      m_z[i]=error_bound[i];
    }
  }
  /** @brief get measurement error bound **/
  std::vector<double> get_measruement_error_bound() {
    std::vector<double> z;
    for(int i=0; i<m_state_alphabet.get_dim(); i++) {
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
#endif /* SYMBOLICMODELGB_HH_ */
