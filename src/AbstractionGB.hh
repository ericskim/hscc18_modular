/*
 * RungeKutta4.hh
 *
 *  created on: 22.04.2015
 *      author: rungger
 */

/** @file **/
#ifndef ABSTRACTIONGB_HH_
#define ABSTRACTIONGB_HH_

#include <iostream>
#include <cstring>
#include <memory>

#include "UniformGrid.hh"
#include "TransitionFunction.hh"
#include "TicToc.hh"

/** @namespace scots **/ 
namespace scots {
/**
 * @class AbstractionGB
 * 
 * @brief Computation of the transition function of symbolic model based on a growth bound
 *
 * See 
 * - the manual in <a href="./../../manual/manual.pdf">manual.pdf</a>
 * - http://arxiv.org/abs/1503.03715 for theoretical background 
 *
 **/
template<class state_type, class input_type>
class AbstractionGB {
public:
  /* @cond  EXCLUDE from doxygen*/
  ~AbstractionGB();

  /* deactivate standard constructor */
  AbstractionGB() = delete;
  /* cannot be copied or moved */
  AbstractionGB(AbstractionGB&&) = delete;
  AbstractionGB(const AbstractionGB&) = delete;
  AbstractionGB& operator=(AbstractionGB&&)=delete;
  AbstractionGB& operator=(const AbstractionGB&)=delete;
  /* @endcond */

  /** @brief constructor with the abstract state alphabet and abstract input
   * alphabet (measurement error bound is set to zero)
   *
   *  @param state_alphabet - UniformGrid representing the abstract state alphabet
   *  @param input_alphabet - UniformGrid representing the abstract input alphabet
   **/
  AbstractionGB(const UniformGrid& state_alphabet,
                const UniformGrid& input_alphabet);

  /** 
   * @brief computes the transition function
   *
   * @param transition_function - the result of the computation
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
   **/
  template<class F1, class F2>
  void compute(TransitionFunction& transition_function,
               F1& system_post,
               F2& radius_post);

  /** 
   * @brief computes the transition function
   *
   * @param transition_function - the result of the computation
   *
   * @param system_post - lambda function as above
   *
   * @param radius_post - lambda function as above
   *
   * @param overflow    - lambda function with signature
   *                      \verbatim [] (abs_type &i) -> bool \endverbatim
   *                      returns true if the abstract state i is in the avoid
   *                      set; otherwise returns false
   *
   **/
  template<class F1, class F2, class F3>
  void compute(TransitionFunction& transition_function,
                    F1& system_post,
                    F2& radius_post,
                    F3&& overflow);
 
  /** @brief get the indices of cells that are used to over-approximated the
   *  attainable set associated with cell (with center x) and input u
   *
   *  @returns std::vector<abs_type>
   *  @param system_post - lambda function as in compute
   *  @param radius_post - lambda function as in compute
   *  @param x - center of cell 
   *  @param u - input
   **/
  template<class F1, class F2>
  std::vector<abs_type> getPost(F1& system_post,
																F2& radius_post,
                                const state_type& x,
                                const input_type& u);


  /** @brief get the center of cells that are used to over-approximated the
   *  attainable set associated with cell with center x and input u
   *
   *  @returns std::vector<state_type>
   *  @param transition_function - the transition function of the abstraction
   *  @param x - center of cell 
   *  @param u - input
   **/
  std::vector<state_type> getPost(const TransitionFunction& transition_function,
                                  const state_type& x,
                                  const input_type& u);

  /** @brief set measurement error bound **/
  void  setMeasurementErrorBound(const state_type& error_bound);	
  /** @brief get measurement error bound **/
  std::vector<double> getMeasruementErrorBound();	
  /** @brief activate console output **/
  void  setVerboseOn();	
  /** @brief deactivate console output **/
  void  setVerboseOff();	
private:
	/** @brief measurement error bound **/
	double* m_z;
  const UniformGrid& m_state_alphabet;
  const UniformGrid& m_input_alphabet;
  /* print progress to the console (default m_verbose=true) */
  bool m_verbose=true;
};

template<class state_type, class input_type>
AbstractionGB<state_type,input_type>::~AbstractionGB() {
  delete[] m_z;
}

template<class state_type, class input_type>
AbstractionGB<state_type,input_type>::AbstractionGB(const UniformGrid& state_alphabet,
                             const UniformGrid& input_alphabet) :
                             m_state_alphabet(state_alphabet),
                             m_input_alphabet(input_alphabet) {
  m_z = new double[state_alphabet.getDimension()];
  for(int i=0; i<state_alphabet.getDimension(); i++) {
    m_z[i]=0;
  }
}

template<class state_type, class input_type>
template<class F1, class F2>
void AbstractionGB<state_type,input_type>::compute(TransitionFunction& transition_function, F1& system_post, F2& radius_post) {
  compute(transition_function,system_post, radius_post, [](const abs_type&) noexcept {return false;});
}

template<class state_type, class input_type>
template<class F1, class F2, class F3>
void AbstractionGB<state_type,input_type>::compute(TransitionFunction& transition_function, F1& system_post, F2& radius_post, F3&& overflow) {
  /* number of cells (=grid points) */
  abs_type N=m_state_alphabet.getTotalNoGridPoints(); 
  /* number of inputs */
  abs_type M=m_input_alphabet.getTotalNoGridPoints();
  /* number of transitions (to be computed) */
  std::size_t T=0; 
  /* state space dimension */
  int dim=m_state_alphabet.getDimension();
  /* for display purpose */
  abs_type counter=0;
  /* some grid information */
  std::vector<abs_type> NN=m_state_alphabet.getNN();
  /* variables for managing the post */
  std::vector<abs_type> lb(dim);  /* lower-left corner */
  std::vector<abs_type> ub(dim);  /* upper-right corner */
  std::vector<abs_type> no(dim);  /* number of cells per dim */
  std::vector<abs_type> cc(dim);  /* coordinate of current cell in the post */
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
    eta[i]=m_state_alphabet.getEta()[i];
    lower_left[i]=m_state_alphabet.getLowerLeftGridPoint()[i];
    upper_right[i]=m_state_alphabet.getUpperRightGridPoint()[i];
  }
  /* number of pre indices of (i,j), initialized to 0  */
  abs_type* no_pre = new abs_type[N*M]();  
  /* number of post indices of (i,j), initialized to 0 */
  abs_type* no_post = new abs_type[N*M]();  
  /* index to pre, where the cell IDs of the pre are stored */
  std::size_t* pre_ptr = new std::size_t[N*M];
  /* lower-left & upper-right corners of hyper rectangle of cells that cover attainable set */
  abs_type* corner_IDs = new abs_type[N*M*2];
  /* is post of (i,j) out of domain ? */
  bool* out_of_domain = new bool[N*M];
  /*
   * first loop: compute corner_IDs:
   * corner_IDs[i*M+j][0] = lower-left cell index of over-approximation of attainable set 
   * corner_IDs[i*M+j][1] = upper-right cell index of over-approximation of attainable set 
   */
  /* loop over all cells */
  for(abs_type i=0; i<N; i++) {
    /* loop over all inputs */
    for(std::size_t j=0; j<M; j++) {
      out_of_domain[i*M+j]=false;
      /* is i an element of the overflow symbols ? */
      if(!j & overflow(i)) {
        for(std::size_t j=0; j<M; j++)
          out_of_domain[i*M+j]=true;
        break;
      }
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
       * covers attainable set 
       */
      abs_type npost=1;
      for(int k=0; k<dim; k++) {
        /* check for out of bounds */
        double left = x[k]-r[k]-m_z[k];
        double right = x[k]+r[k]+m_z[k];
        if(left <= lower_left[k]-eta[k]/2.0  || right >= upper_right[k]+eta[k]/2.0) {
          out_of_domain[i*M+j]=true;
          break;
        } 
        /* integer coordinate of lower left corner of post */
        lb[k] = static_cast<abs_type>((left-lower_left[k]+eta[k]/2.0)/eta[k]);
        /* integer coordinate of upper right corner of post */
        ub[k] = static_cast<abs_type>((right-lower_left[k]+eta[k]/2.0)/eta[k]);
        /* number of grid points in the post in each dimension */
        no[k]=(ub[k]-lb[k]+1);
        /* total number of post */
        npost*=no[k];
        cc[k]=0;
      }
      corner_IDs[i*(2*M)+2*j]=0;
      corner_IDs[i*(2*M)+2*j+1]=0;
      if(out_of_domain[i*M+j]) 
        continue;
      /* compute indices of post */
      for(abs_type k=0; k<npost; k++) {
        abs_type q=0;
        for(int l=0; l<dim; l++) 
          q+=(lb[l]+cc[l])*NN[l];
        cc[0]++;
        for(int l=0; l<dim-1; l++) {
          if(cc[l]==no[l]) {
            cc[l]=0;
            cc[l+1]++;
          }
        }
        /* (i,j,q) is a transition */    
        /* increment number of pres for (q,j) */ 
        no_pre[q*M+j]++;
        /* store id's of lower-left and upper-right cell */
        if(k==0)
          corner_IDs[i*(2*M)+2*j]=q;
        if(k==npost-1)
          corner_IDs[i*(2*M)+2*j+1]=q;
      }
      /* increment number of transitions by number of post */
      T+=npost;
      no_post[i*M+j]=npost;
    }
    /* print progress */
    if(m_verbose && ((double)i/(double)N*100)>counter){
      if(counter==0)
        std::cout << "1st loop: ";
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
  /* compute pre_ptr */
  abs_type sum=0;
  for(abs_type i=0; i<N; i++) {
    for(abs_type j=0; j<M; j++) {
      sum+=no_pre[i*M+j];
      pre_ptr[i*M+j]=sum;
    }
  }
  /* allocate memory for pre list */
  abs_type* pre = new abs_type[T];

  /*
   * second loop: fill pre array
   */
  counter=0;
  for(abs_type i=0; i<N; i++) {
    /* loop over all inputs */
    for(abs_type j=0; j<M; j++) {
    /* is x an element of the overflow symbols ? */
      if(out_of_domain[i*M+j]) 
        continue;
      /* extract lower-left and upper-bound points */
      abs_type k_lb=corner_IDs[i*2*M+2*j];
      abs_type k_ub=corner_IDs[i*2*M+2*j+1];
      abs_type npost=1;
      /* cell idx to coordinates */
      for(int k=dim-1; k>=0; k--) {
        /* integer coordinate of lower left corner */
        lb[k]=k_lb/NN[k];
        k_lb=k_lb-lb[k]*NN[k];
        /* integer coordinate of upper right corner */
        ub[k]=k_ub/NN[k];
        k_ub=k_ub-ub[k]*NN[k];
        /* number of grid points in each dimension in the post */
        no[k]=(ub[k]-lb[k]+1);
        /* totoal no of post of (i,j) */
        npost*=no[k];
        cc[k]=0;
      }
      for(abs_type k=0; k<npost; k++) {
        abs_type q=0;
        for(int l=0; l<dim; l++) 
          q+=(lb[l]+cc[l])*NN[l];
        cc[0]++;
        for(int l=0; l<dim-1; l++) {
          if(cc[l]==no[l]) {
            cc[l]=0;
            cc[l+1]++;
          }
        }
        /* (i,j,q) is a transition */
        pre[--pre_ptr[q*M+j]]=i;
      }
    }
    /* print progress */
    if(m_verbose && ((double)i/(double)N*100)>counter){
      if(counter==0)
        std::cout << "2nd loop: ";
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
  /* set values of abstract system */
  transition_function.m_no_states=N;
  transition_function.m_no_inputs=M;
  transition_function.m_no_transitions=T;
  transition_function.m_no_pre=no_pre;
  transition_function.m_no_post=no_post;
  transition_function.m_pre=pre;
  transition_function.m_pre_ptr=pre_ptr;
  /* cleanup */
  delete[] corner_IDs;
  delete[]  out_of_domain;
}

template<class state_type, class input_type>
std::vector<state_type>  AbstractionGB<state_type,input_type>::getPost(const TransitionFunction& abstraction, const state_type& x, const input_type& u) {
  abs_type i,j;
  std::vector<state_type> post;
  std::vector<abs_type> k;

  m_state_alphabet.xtoi(i,x);
  m_input_alphabet.xtoi(j,u);
  k=abstraction.getPost(i,j);
  for(abs_type v=0; v<k.size(); v++) {
    state_type s;
    m_state_alphabet.itox(k[v],s);
    post.push_back(s);
  }
  return std::move(post);
}

template<class state_type, class input_type>
template<class F1, class F2> 
std::vector<abs_type> AbstractionGB<state_type,input_type>::getPost(F1& system_post, F2& radius_post, const state_type& x, const input_type& u) {
	/* initialize return value */
  std::vector<abs_type> post;
  post.clear();
  /* state space dimension */
	int dim = m_state_alphabet.getDimension();
  /* variables for managing the post */
  std::vector<abs_type> lb(dim);  /* lower-left corner */
  std::vector<abs_type> ub(dim);  /* upper-right corner */
  std::vector<abs_type> no(dim);  /* number of cells per dim */
  std::vector<abs_type> cc(dim);  /* coordinate of current cell in the post */
  /* get radius */
	state_type r;
  state_type eta;
  /* for out of bounds check */
  state_type lower_left;
  state_type upper_right;
  /* fill in data */
	for(int i=0; i<dim; i++) {
    eta[i]=m_state_alphabet.getEta()[i];
    r[i]=m_state_alphabet.getEta()[i]/2.0+m_z[i];
    lower_left[i]=m_state_alphabet.getLowerLeftGridPoint()[i];
    upper_right[i]=m_state_alphabet.getUpperRightGridPoint()[i];
  }

  /* make a copy of cell */
	state_type xx=x;

  /* compute growth bound and numerical solution of ODE */
  radius_post(r,xx,u);
  system_post(x,u);

  /* determine post */
  abs_type npost=1;
  for(int k=0; k<dim; k++) {
    /* check for out of bounds */
    double left = x[k]-r[k]-m_z[k];
    double right = x[k]+r[k]+m_z[k];
    if(left <= lower_left[k]-eta[k]/2.0  || right >= upper_right[k]+eta[k]/2.0) {
      return post;
    } 
    /* integer coordinate of lower left corner of post */
    lb[k] = static_cast<abs_type>((left-lower_left[k]+eta[k]/2.0)/eta[k]);
    /* integer coordinate of upper right corner of post */
    ub[k] = static_cast<abs_type>((right-lower_left[k]+eta[k]/2.0)/eta[k]);
    /* number of grid points in the post in each dimension */
    no[k]=(ub[k]-lb[k]+1);
    /* total number of post */
    npost*=no[k];
    cc[k]=0;
  }
  /* compute indices of post */
  for(abs_type k=0; k<npost; k++) {
    abs_type q=0;
    for(int l=0; l<dim; l++)  {
      q+=(lb[l]+cc[l])*m_state_alphabet.getNN()[l];
    }
    cc[0]++;
    for(int l=0; l<dim-1; l++) {
      if(cc[l]==no[l]) {
        cc[l]=0;
        cc[l+1]++;
      }
    }
    /* (i,j,q) is a transition */    
    post.push_back(q);
  }
  return post;
}

template<class state_type, class input_type>
void  AbstractionGB<state_type,input_type>::setMeasurementErrorBound(const state_type& error_bound) {
  for(int i=0; i<m_state_alphabet.getDimension(); i++) {
		m_z[i]=error_bound[i];
  }
}

template<class state_type, class input_type>
std::vector<double> AbstractionGB<state_type,input_type>::getMeasruementErrorBound() {
  std::vector<double> z;
  for(int i=0; i<m_state_alphabet.getDimension(); i++) {
		z.push_back(m_z[i]);
  }
	return z;
}

template<class state_type, class input_type>
void AbstractionGB<state_type,input_type>::setVerboseOn() {
  m_verbose=true;
}

template<class state_type, class input_type>
void AbstractionGB<state_type,input_type>::setVerboseOff() {
  m_verbose=false;
}

} /* close namespace */

#endif /* ABSTRACTIONGB_HH_ */