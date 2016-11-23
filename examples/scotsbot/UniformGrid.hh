/*
 * UniformGrid.hh
 *
 *  created on: 26.10.2015
 *      author: rungger
 */

#ifndef UNIFORMGRID_HH_
#define UNIFORMGRID_HH_

#include <vector>
#include <sstream>
#include <string>
#include <cassert>
#include <cmath>
#include <set>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdexcept>

namespace scots {

/* class: UnifromGrid
 *
 * stores information of a uniform grid
 *
 *
 * Properties: 
 * - the grid points are distributied uniformly in each dimension 
 * - the domain of the unfiorm grid is defined by a hyper interval
 * - each grid pont is associated with a cell, i.e. a hyper rectangle with
 *   radius eta/2+z centered at the grid point 
 * 
 * Grid point alignment:
 * - the origin is a grid point (not necessarily contained in the set)
 * - the distance of the grid points in each dimension i is defined by eta[i]
 *
 * See 
 * - http://arxiv.org/abs/1503.03715 for theoretical background
 *
 */
template<class grid_point_t>
class UniformGrid {
  friend class ReachabilityGame;
  friend class IO;
  private:
    /* var: dim_ 
     * dimension of the real space */
    size_t dim_;
    /* var: eta_
     * dim_-dimensional vector containing the grid node distances */
    grid_point_t eta_;
    /* var: z_
     * dim_-dimensional vector containing the measurement error bound */
    grid_point_t z_;
    /* var: firstGridPoint_
     * dim_-dimensinal vector containing the real values of the first grid point */
    grid_point_t firstGridPoint_;
    /* var: nofGridPoints_
     * integer array[dim_] containing the number of grid points in each dimension */
    std::vector<size_t> nofGridPoints_;
    /* var: N_ 
     * total number of grid points */
    size_t N_;
    std::vector<size_t> NN_;

    /* var: abstractSet_ 
     * contains a list of indices with are in the set */
    std::set<size_t> abstractSet_;

  public:
    /* constructor: UniformGrid 
     * provide uniform grid parameters and domain defining hyper interval
     * 
     * Input:
     * lb    - lower left corner of the domain hyper interval
     * ub    - upper right corner of the domain hyper interval
     * eta   - grid point distances
     *
     */
    UniformGrid(const size_t dim, const grid_point_t &lb, const grid_point_t &ub, const grid_point_t &eta) {
      for (size_t i=0; i<dim; i++) {
        if((lb[i] > ub[i])) {
            std::ostringstream os;
            os << "Error: scots::UniformGrid: lower bound must be less than or equal upper bound.";
            throw std::invalid_argument(os.str().c_str());
        }
      } 
      /* measurement error bound is set to zero */
      grid_point_t z;
      for(size_t i=0; i<dim; i++)
        z[i]=0;
      initGrid(dim,lb,ub,eta,z);
    }
    /* constructor: UniformGrid 
     * provide uniform grid parameters and domain defining hyper interval
     * 
     * Input:
     * lb    - lower left corner of the domain hyper interval
     * ub    - upper right corner of the domain hyper interval
     * eta   - grid point distances
     * z     - measurement error bound 
     */
    UniformGrid(const size_t dim, const grid_point_t &lb, const grid_point_t &ub, const grid_point_t &eta, const grid_point_t &z) {
      for (size_t i=0; i<dim; i++) {
        if((lb[i] > ub[i])) {
            std::ostringstream os;
            os << "Error: scots::UniformGrid: lower bound must be less than or equal upper bound.";
            throw std::invalid_argument(os.str().c_str());
        }
      }
      initGrid(dim,lb,ub,eta,z);
    }

    /* constructor: UniformGrid
     * empty
     */
    UniformGrid() {
    }

    /* constructor: UniformGrid
     * 
     * copy constructor
     */
    UniformGrid(const UniformGrid& other) {
     *this=other;
    }

    ~UniformGrid() {
    }
   
    /* function: copy assignment operator */
    UniformGrid& operator=(const UniformGrid &other) {
      dim_=other.dim_;
      N_=other.N_;
      z_ = other.z_;
      eta_ = other.eta_;
      firstGridPoint_ = other.firstGridPoint_;
      NN_ = other.NN_;
      return *this;
    }
    /* function: xtoi
     * compute the index idx associated with a  */
    inline void xtoi(size_t &idx, const grid_point_t x) const {
      long int t;
      idx=0;
      for(size_t k=0; k<dim_; k++) {
       t = lround(((x[k]-firstGridPoint_[k])/eta_[k]));
       if (t< 0 || t > (long int)nofGridPoints_[k]) {
          std::ostringstream os;
          os << "Error: UniformGrid: x is outside uniform grid.";
          throw std::runtime_error(os.str().c_str());
       }
       idx+=t*NN_[k];
      }
    }
    /* function: itox
    * compute the cooridnated of the grid point associated with the index idx */
    inline void itox(size_t idx, grid_point_t &x) const {
      if (idx > N_){
          std::ostringstream os;
          os << "Error: UniformGrid: idx larger than number of grid points.";
          throw std::runtime_error(os.str().c_str());
      }
      size_t num;
      for(int i=dim_-1; i>=0; i--) {
        if(i) {
          num=idx/NN_[i];
          idx=idx%NN_[i];
        } else {
          num=idx;
        }
        x[i]=firstGridPoint_[i]+num*eta_[i];
      }
    }
    /* function: printInfo
     * print some numbers related to the symbolic set*/
    void printInfo(int verbosity=0) const {
      std::cout << "Grid node distance (eta)in each dimension: ";
      for(size_t i=0; i<dim_; i++) 
        std::cout << eta_[i] << " ";
      std::cout << std::endl;
      if (verbosity) {
        std::cout << "First grid point: ";
        for(size_t i=0; i<dim_; i++) 
            std::cout << firstGridPoint_[i] << " ";
          std::cout << std::endl;
      }
      std::cout << "Number of grid points in each dimension: ";
      for(size_t i=0; i<dim_; i++) 
        std::cout << nofGridPoints_[i] << " ";
      std::cout << std::endl;

      std::cout << "Number of grid points: "<< N_ << std::endl;
    }

    /* function: remIndices
     *
     * remove all indices i from the abstract set
     * for which the function set(i) == true 
     *
     */
    template<class F>
    void remIndices(F &set) {
      /* iterate over elements of the set */
      std::set<size_t>::iterator it;
      for (it=abstractSet_.begin(); it != abstractSet_.end(); ) {
        if(set(*it))
          it = abstractSet_.erase(it);
        else
          ++it;
      }
    } 
    /* function: addIndices
     *
     * add all indices i to the abstract set for which the function 
     * set(i) == true 
     */
    template<class F>
    void addIndices(F &set) {
      for(size_t i=0; i<N_; i++) {
        if(set(i))
          abstractSet_.insert(i);
      }
    }
    /* function: remGridPoints
     *
     * remove all grid points x from the abstract set
     * for which the function set(x) == true 
     *
     */
    template<class F>
    void remGridPoints(F &set) {
      grid_point_t x;
      /* iterate over elements of the set */
      std::set<size_t>::iterator it;
      for (it=abstractSet_.begin(); it != abstractSet_.end(); ) {
        itox(*it,x);
        if(set(x))
          it = abstractSet_.erase(it);
        else
          ++it;
      }
    } 
    /* function: addGridPoints
     *
     * add all grid points x to the abstract set for which the function 
     * set(x) == true 
     */
    template<class F>
    void addGridPoints(F &set) {
      grid_point_t x;
      for(size_t i=0; i<N_; i++) {
        itox(i,x);
        if(set(x))
          abstractSet_.insert(i);
      }
    }
    /* function: fillAbstractSet 
     * add all grid points to the abstract set
     */
    void fillAbstractSet() {
      for(size_t i=0; i<N_; i++) 
        abstractSet_.insert(i);
    }
    /* function: clearAbstractSet 
     * remove all grid points to the abstract set
     */
    void clearAbstractSet() {
      for(size_t i=0; i<N_; i++) 
        abstractSet_.clear();
    }
    /* function:  getDimension 
     * get the dimension of the real space of the uniform grid */
    inline size_t getDimension(void) const {
      return dim_;
    }
    /* function:  getZ
     * returns z_*/
    inline const grid_point_t getZ() const {
      return z_;
    }
    /* function:  getEta
     * return  eta_*/
    inline const grid_point_t getEta() const {
      return eta_;
    }
    /* function:  getFirstGridPoint 
     * returns the first grid point */
    inline const grid_point_t getFirstGridPoint() const {
      return firstGridPoint_;
    }
    /* function:  getNofGridPoints
     * return pointer to size_t array containing the number of grid points*/
    inline const std::vector<size_t> getNofGridPoints() const {
        return nofGridPoints_;
    }  
    /* function:  getSetSize 
     * get the number of elements in the abstract set */
    inline size_t getSetSize(void) const {
      return abstractSet_.size();
    }
    /* function:  getN
     * return  number of grid points*/
    inline size_t getN() const {
        return N_;
    }
    inline const std::vector<size_t> getNN() const {
        return NN_;
    }
    /* function: setIdxToDoubleMatrix
     * the grid points in the abstract set are saved to the double array of size (dim_) x (setSize) */
    void setIdxToDoubleMatrix(double gridPoints[]) const {
      size_t num;
      size_t k=0;
      size_t K=abstractSet_.size();
      std::set<size_t>::iterator it;
      for (it=abstractSet_.begin(); it != abstractSet_.end();++it ) {
        size_t idx = *it;
        for(int i=dim_-1; i>=0; i--) {
          num=idx/NN_[i];
          idx=idx%NN_[i];
          gridPoints[i*K+k]=firstGridPoint_[i]+num*eta_[i];
        }
        k++;
      }
    }

    /* function: projectSet
     * the grid points in abstract set are projected onto the specified dimensions
     * projectDimension: 1,2,3 ... dim_ */
    std::vector<size_t> projectSet(std::vector<size_t> projectDimension) {
      /* check input data */
      size_t n=projectDimension.size();
      if(n>dim_) {
        std::ostringstream os;
        os << "Error: UniformGrid: project dimension larger than dimension of grid points.";
        throw std::runtime_error(os.str().c_str());
      }
      for(size_t i=0;i<n;i++) {
        if(projectDimension[i]>dim_) {
          std::ostringstream os;
          os << "Error: UniformGrid: Cannot project UniformGrid onto given dimensions.";
          throw std::runtime_error(os.str().c_str());
        }
      }

      std::vector<size_t> result;
      result.clear();

      size_t num;
      size_t cor[dim_];

      /* loop over all the elements in abstractSet */
      std::set<size_t>::iterator it;
      for(it=abstractSet_.begin(); it != abstractSet_.end();++it) {
        size_t idx=*it;
        /* find out the coordinate of each element */
        for(int j=dim_-1; j>=0; j--) {
          num=idx/NN_[j];
          idx=idx%NN_[j];
          cor[j]=num;
        }
        /* find out the projected coordinate of each element */
        size_t* corProj=(size_t*)calloc(dim_,sizeof(size_t));
        for(size_t k=0;k<n;k++) {
          corProj[projectDimension[k]-1]=cor[projectDimension[k]-1];
        }
        /* calculate the projected idx */
        size_t idxProj=0;
        for(size_t v=0;v<dim_;v++) {
          idxProj+=corProj[v]*NN_[v];
        }
        free(corProj);
        /* check if the projected idx is already in the result set */
        if(!(std::find(result.begin(),result.end(),idxProj)!=result.end()))
          result.push_back(idxProj);
      }

      return result;

    }

    /* function: project
     * the grid points specified by argument idxSet are projected onto the specified dimensions
     * projectDimension: 1,2,3 ... dim_ */
    std::vector<size_t> project(std::vector<size_t> idxSet, std::vector<size_t> projectDimension) {
      /* check input data */
      size_t n=projectDimension.size();
      if(n>dim_) {
        std::ostringstream os;
        os << "Error: UniformGrid: project dimension larger than dimension of grid points.";
        throw std::runtime_error(os.str().c_str());
      }
      for(size_t i=0;i<n;i++) {
        if(projectDimension[i]>dim_) {
          std::ostringstream os;
          os << "Error: UniformGrid: Cannot project UniformGrid onto given dimensions.";
          throw std::runtime_error(os.str().c_str());
        }
      }

      std::vector<size_t> result;
      result.clear();

      size_t num;
      size_t cor[dim_];

      /* loop over all the elements in idxSet */
      for(size_t i=0;i<idxSet.size();i++) {
        size_t idx=idxSet[i];
        /* find out the coordinate of each element */
        for(int j=dim_-1; j>=0; j--) {
          num=idx/NN_[j];
          idx=idx%NN_[j];
          cor[j]=num;
        }
        /* find out the projected coordinate of each element */
        size_t* corProj=(size_t*)calloc(dim_,sizeof(size_t));
        for(size_t k=0;k<n;k++) {
          corProj[projectDimension[k]-1]=cor[projectDimension[k]-1];
        }
        /* calculate the projected idx */
        size_t idxProj=0;
        for(size_t v=0;v<dim_;v++) {
          idxProj+=corProj[v]*NN_[v];
        }
        free(corProj);
        /* check if the projected idx is already in the result set */
        if(!(std::find(result.begin(),result.end(),idxProj)!=result.end()))
          result.push_back(idxProj);
      }

      return result;

    }

  private:
    /* determine the number of grid points and fix first/last grid point */
    void initGrid(const size_t dim, const grid_point_t &lb, const grid_point_t &ub, const grid_point_t &eta, const grid_point_t &z) {
      dim_=dim;
      /* init integer arrays */
      nofGridPoints_.resize(dim);
      NN_.resize(dim);
      /* determine number grid points in each dimension */
      double Nl, Nu;
      for (size_t i=0; i<dim; i++) {
        z_[i]=z[i];
        eta_[i]=eta[i];

        Nl=std::ceil(lb[i]/eta[i]);
        Nu=std::floor(ub[i]/eta[i]);
        /* number of grid points */
        nofGridPoints_[i]=(size_t)std::abs(Nu-Nl)+1;
        /* first grid point coordinates */
        firstGridPoint_[i]=Nl*eta[i];
      }
      /* total number of total grid points */
      N_=1;
      for(size_t i=0; i<dim_; i++) {
          NN_[i]=N_;
          N_*=nofGridPoints_[i];
      }
    }
}; /* close class def */

} /* close namespace */


#endif /* UNIFORMGRID_HH_ */ 
