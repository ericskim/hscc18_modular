/*
 * IndexSet.hh
 *
 *  created: Jan 2016
 *   author: Matthias Rungger
 */

/** @file **/

#ifndef INDEXSET_HH_
#define INDEXSET_HH_

#include <vector>
#include <iostream>

/* cudd library */
#include "cuddObj.hh"

namespace scots {

/**
 * @brief id_type (default = uint32_t) is used to represent the indices in the IndexSet\n
 *        it should be equal abs_type if IndexSet is used in connection with
 *        UniformGrid
 **/
using id_type = std::uint32_t;

/** @class IndexSet 
 *
 *  @brief A BDD representation of a multidimensional index set
 * 
 * The index set with dimension \f$d \in \mathbb N \f$ is defined by the number of
 * indices \f$i_0,\ldots,i_{d-1}\f$ and given by
 *
 * \f$ S:= \{0,\ldots,i_0-1\} \times  \cdots  \times \{0, \ldots, i_{d-1}-1\}\f$
 *
 * each index vector \f$ s \in S \f$ is associated with an ID 
 * \f$ \mathrm{ID} = \prod_{j=0}^{d-1}i_j s_j \f$
 *
 **/
class IndexSet {
private:
  /* reference to BDD manager */
  const Cudd& m_manager;
  /* dimension of the index set */
  int m_dim=0;
  /* m_size = i_0 * .. * i_{m_dim-1} */
  id_type m_size=0;
  /* no of indices in each dimension */
  std::vector<id_type> m_no_indices;
  /* each vector m_bdd_var[d] contains the BDD variables dimension d */
  std::vector<std::vector<BDD>> m_bdd_var;
  /* each vector m_bdd_var_id[d] contains the BDD variable IDs in dimension d */
  std::vector<std::vector<unsigned int>> m_bdd_var_id;
  /* an array[m_size] containing the BDD representation of each element */
  std::vector<std::vector<BDD>> m_id_to_bdd;
  /* total number of bdd vars */
  unsigned int m_no_bdd_var;
  /* helpder m_NN */
  std::vector<id_type> m_NN;
  /* helper function */
  void add_one(int* phase, unsigned int num) {
    int carry = 1;
    for(int i=num-1; carry && (i>=0) ;i--) {
      if(phase[i]) {
        phase[i]=0;
        carry=1;
      } else {
        phase[i]=1;
        carry=0;
      }
    }
  }

public:
  /** @cond  EXCLUDE from doxygen **/
  /* default destructor */
  ~IndexSet()=default;
  /* deactivate move constructor */
  IndexSet(IndexSet&&)=delete;
  /* deactivate move assignment operator */
  IndexSet& operator=(IndexSet&&)=delete;
  /* copy constructor deleted (cannot be copied) */
  IndexSet(const IndexSet&)=delete;
  /* copy assignment operator */
  IndexSet& operator=(const IndexSet&)=delete;
  /* @endcond */
  
  IndexSet(const Cudd& manager,
           std::vector<id_type> no_indices) :
           m_manager(manager), m_no_indices(no_indices) {
    m_dim=m_no_indices.size();
    /* compute size */ 
    m_size=1;
    m_NN.resize(m_dim);
    m_NN[0]=1;
    for(int i=0; i<m_dim; i++) {
      /* check overflow */
      if(m_size > (std::numeric_limits<id_type>::max()/m_no_indices[i])) {
        throw std::runtime_error("\nscots::IndexSet: number of indices exceeds maximum value of id_type.");
      }
      if(i) {
        m_NN[i]=m_no_indices[i-1]*m_NN[i-1];
      }
      m_size*=m_no_indices[i];
    }
    /* compute no of BDD neccessary to represent the set index set per dim */
    m_no_bdd_var = 0;
    std::vector<int> no_bdd_var_per_dim(m_dim,0);
    for(int i=0; i<m_dim; i++) {
      id_type x = no_indices[i]-1u;
      while(x) {
        x>>=1u;
        no_bdd_var_per_dim[i]++;
      }
      m_no_bdd_var+=no_bdd_var_per_dim[i];
    }
    /* create new BDD variables to represent the IndexSet */
    m_bdd_var.resize(m_dim);
    m_bdd_var_id.resize(m_dim);
    /* compute m_id_to_bdd array for each coordinate */
    m_id_to_bdd.resize(m_dim);
    for(int i=0; i<m_dim; i++) {
      int* phase = new int[no_bdd_var_per_dim[i]] (); 
      BDD* vars = new BDD[no_bdd_var_per_dim[i]];
      /* get new BDD variable IDs */
      for(int j=0; j<no_bdd_var_per_dim[i]; j++) {
        vars[j] = manager.bddVar();
        m_bdd_var[i].push_back(vars[j]);
        m_bdd_var_id[i].push_back(vars[j].NodeReadIndex());
      }
      m_id_to_bdd[i].resize(m_no_indices[i]);
      for(id_type id=0; id<m_no_indices[i]; id++) {
        m_id_to_bdd[i][id] = manager.bddComputeCube(vars,phase,no_bdd_var_per_dim[i]);
        add_one(phase,no_bdd_var_per_dim[i]);
        //m_id_to_bdd[i][id].PrintMinterm();
      }
      delete[] vars;
      delete[] phase;
    }
  }

  BDD id_to_bdd(id_type id) const {
    BDD bdd = m_manager.bddOne();
    id_type num;
    for(int k=m_dim-1; k > 0; k--) {
      num=id/m_NN[k];
      id=id%m_NN[k];
      bdd = bdd & m_id_to_bdd[k][num];
    }
    num=id;
    bdd = bdd & m_id_to_bdd[0][num];
    return bdd;
  }

  BDD interval_to_bdd(const std::vector<id_type>& lb,   
                      const std::vector<id_type>& ub) const {
    BDD bdd = m_manager.bddOne();
    for(int i=0; i<m_dim; i++) {
      unsigned int lowerB=static_cast<unsigned int>(lb[i]);
      unsigned int upperB=static_cast<unsigned int>(ub[i]);
      bdd = bdd & m_manager.Interval(m_bdd_var[i], lowerB, upperB);
    }
    return bdd;
  }


  void print_info(int verbose=0) {
    std::cout << "Number of indices: "<< m_size << "\n";
    std::cout << "Number of indices per dimension ";
    for(int i=0; i<m_dim; i++) {
      std::cout << m_no_indices[i] << " ";
    }
    std::cout << "\nNumber of BDD variables per dimension ";
    for(int i=0; i<m_dim; i++) {
      std::cout << m_bdd_var_id[i].size() << " ";
    }
    if(verbose) {
      std::cout << "\nBDD variable IDs: ";
      for(int i=0; i<m_dim; i++) {
        std::cout << "\n in dimension " << i << " : ";
        for(size_t j=0; j<m_bdd_var_id[i].size(); j++) 
        std::cout << m_bdd_var_id[i][j] << " ";
      }
    }
    std::cout << "\n\n";
  }


 
//
//
//    int storeReturnValue = Dddmp_cuddBddStore(
//      mdest.getManager(),
//      NULL,
//      tosave.getNode(),
//      //(char**)varnameschar, // char ** varnames, IN: array of variable names (or NULL)
//      NULL, // char ** varnames, IN: array of variable names (or NULL)
//      NULL,
//      DDDMP_MODE_BINARY,
//      // DDDMP_VARNAMES,
//      DDDMP_VARIDS,
//      NULL,
//      file
//    );
//
//    fclose(file);
//    if (storeReturnValue!=DDDMP_SUCCESS) 
//      throw "Error: Unable to write BDD to file.";
//    else
//      std::cout << "Symbolic set saved to file: "<< filename << std::endl;

  BDD get_zero() const {
    return m_manager.bddZero();
  }

  BDD get_one() const {
    return m_manager.bddOne();
  }

    return m_manager.getManager();
  }



  unsigned int get_no_bdd_var(void) const {
    return m_no_bdd_var;
  }

}; /* close class def */
} /* close namespace */
#endif /* INDEXSET_HH_ */
