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
public:
  /* reference to BDD manager */
  const Cudd& m_manager;
  /* m_size = i_0 * .. * i_{m_dim-1} */
  id_type m_size=0;
  /* no of indices in each dimension */
  std::vector<id_type> m_no_indices;
  /* each vector m_bdd_var[d] contains the BDD variable IDs in dimension d */
  std::vector<std::vector<unsigned int>> m_bdd_var;
  /* an array[m_size] containing the BDD representation of each element */
  std::vector<BDD> m_id_to_bdd;
  unsigned int m_no_bdd_var;
  /* helper function */
  void add_one(int* phase) {
    int carry = 1;
    for(unsigned int i=0;carry;i++) {
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
    int dim=m_no_indices.size();
    /* compute size */ 
    m_size=1;
    for(int i=0; i<dim; i++) {
      /* check overflow */
      if(m_size > (std::numeric_limits<id_type>::max()/m_no_indices[i])) {
        throw std::runtime_error("\nscots::IndexSet: number of indices exceeds maximum value of id_type.");
      }
      m_size*=m_no_indices[i];
    }
    /* compute no of BDD neccessary to represent the set index set per dim */
    m_no_bdd_var = 0;
    std::vector<int> no_bdd_var_per_dim(dim,0);
    for(int i=0; i<dim; i++) {
      id_type x = no_indices[i]-1u;
      while(x) {
        x>>=1u;
        no_bdd_var_per_dim[i]++;
      }
      m_no_bdd_var+=no_bdd_var_per_dim[i];
    }
    BDD* vars = new BDD[m_no_bdd_var];
    /* create new BDD variables to represent the IndexSet */
    m_bdd_var.resize(dim);
    unsigned int shift = 0;
    for(int i=0; i<dim; i++) {
      if(i) {
        shift += no_bdd_var_per_dim[i-1];
      }
      for(int j=0; j<no_bdd_var_per_dim[i]; j++) {
        vars[shift+j] = manager.bddVar();
        m_bdd_var[i].push_back(vars[shift+j].NodeReadIndex());
      }
    }
    /* compute m_id_to_bdd array */
    int* phase = new int[m_no_bdd_var] (); 
    m_id_to_bdd.resize(m_size);
    for(id_type id=0; id<m_size; id++) {
      m_id_to_bdd[id] = manager.bddComputeCube(vars,phase,m_no_bdd_var);
      add_one(phase);
    }
    delete[] vars;
    delete[] phase;
  }

  void print_info(int verbose=0) {
    int dim = m_no_indices.size();
    std::cout << "Number of indices: "<< m_size << "\n";
    std::cout << "Number of indices per dimension ";
    for(int i=0; i<dim; i++) {
      std::cout << m_no_indices[i] << " ";
    }
    std::cout << "\nNumber of BDD variables per dimension ";
    for(int i=0; i<dim; i++) {
      std::cout << m_bdd_var[i].size() << " ";
    }
    if(verbose) {
      std::cout << "\nBDD variable IDs: ";
      for(int i=0; i<dim; i++) {
        std::cout << "\n in dimension " << i << " : ";
        for(size_t j=0; j<m_bdd_var[i].size(); j++) 
        std::cout << m_bdd_var[i][j] << " ";
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

}; /* close class def */
} /* close namespace */
#endif /* INDEXSET_HH_ */
