/*
 * UniformGrid.hh
 *
 *     created: Dec 2016
 *      author: Matthias Rungger
 *              Frederik Kunik
 */

/** @file **/

#ifndef UNIFORMGRID_HH_
#define UNIFORMGRID_HH_

#include <vector>
#include <array>
#include <cmath>
#include <sstream>
#include <string>
#include <exception>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <limits>

#include "TransitionFunction.hh"
//#include "FileHandler.hh"

#define SCOTS_UG_TYPE   "UNIFORMGRID"
#define SCOTS_UG_DIM    "DIM"
#define SCOTS_UG_FIRST  "FIRST"
#define SCOTS_UG_ETA    "ETA"
#define SCOTS_UG_GPPD   "GPPD"

/** @namespace scots **/ 
namespace scots {

/**
 * @class UniformGrid 
 *
 * @brief Holds the information of a uniform grid confined by a hyper-interval
 *
 * Properties: 
 * - the grid points are distributied uniformly in each dimension 
 * - the domain of the unfiorm grid is defined by a hyper interval 
 * - each grid pont is associated with a hyper-rectangle with radius eta/2 centered at the grid point 
 *
 * Grid point alignment: 
 * - the origin is a grid point (not necessarily contained in the set) 
 * - the distance of the grid points in each dimension i is defined by eta[i] 
 *
 * See 
 * - the manual in <a href="./../../manual/manual.pdf">manual</a>
 * - http://arxiv.org/abs/1503.03715 for theoretical background 
 **/
class UniformGrid {
private:
  /** @brief dimension of the Eucleadian space **/
  int m_dim;                
  /** @brief m_dim-dimensional vector containing the grid node distances **/
  double* m_eta;                      
	/** @brief m_dim-dimensional vector containing the real values of the first grid point **/
  double* m_first;        
	/** @brief scots::abs_type array[m_dim] containing the number of grid points in each dimension **/
  abs_type* m_no_grid_points;        
	/** @brief array recursively defined by: m_NN[0]=1; m_NN[i]=m_NN[i-1}*no_grid_points[i-1]; **/
  abs_type* m_NN;                       

public:
  /* @cond  EXCLUDE from doxygen */
  /* default constructor */
  UniformGrid() {
    m_dim = 0;
    m_eta = nullptr;
    m_first = nullptr;
    m_no_grid_points = nullptr; 
    m_NN = nullptr;
  }
  /* destructor */
  ~UniformGrid() {
    delete[] m_eta;
    delete[] m_first;
    delete[] m_no_grid_points;
    delete[] m_NN;
  }
  /* copy constructor */
  UniformGrid(const UniformGrid& other) : UniformGrid() {
    *this=other;
  }
  /* move constructor */
  UniformGrid(UniformGrid&& other) : UniformGrid() {
    *this=std::move(other);
  }
  /* copy assignment operator */
  UniformGrid& operator=(const UniformGrid &other) {
    if(this==&other)
      return *this;
    reset();
    m_dim=other.m_dim;
    if(m_dim != 0) {
      m_eta = new double[m_dim];
      m_first = new double[m_dim];
      m_no_grid_points = new abs_type[m_dim];
      m_NN = new abs_type[m_dim];
      for(int index=0; index<m_dim; index++) {
        m_eta[index] = other.m_eta[index];
        m_first[index]  = other.m_first[index];
        m_no_grid_points[index]  = other.m_no_grid_points[index];
        m_NN[index]  = other.m_NN[index];
      }
    } 
    return *this;
  }
  /* move assignment operator */
  UniformGrid& operator=(UniformGrid&& other) {
    reset();

    m_dim=other.m_dim;
    m_eta=other.m_eta;
    m_first=other.m_first;
    m_no_grid_points=other.m_no_grid_points;
    m_NN=other.m_NN;

    other.null();
    return *this;
  } 
  /* @endcond */

  /**
   * @brief provide uniform grid parameters and domain defining hyper interval 
   * 
   * @param dim   - dimension of the real space
   * @param lb    - lower-left corner of the hyper-interval confining the uniform grid
   * @param ub    - upper-right corner of the hyper-interval confining the uniform grid
   * @param eta   - grid point distances
   *
   */
  template<class grid_point_t>
  UniformGrid::UniformGrid(const int dim,
                           const grid_point_t& lb,
                           const grid_point_t& ub,
                           const grid_point_t& eta) : UniformGrid() {
    m_dim = dim;
    if(m_dim != 0) {
      /* check inut arguments */
      for(int index=0; index<dim; index++) {
        if(eta[index] <= 0) 
          throw std::runtime_error("\nscots::UniformGrid: eta must have positive entries.");
        if(lb[index] > ub[index]) 
          throw std::runtime_error("\nscots::UniformGrid: lower-left bound must be less than or equal to upper-right bound.");
      }
      m_eta = new double[m_dim];
      m_first = new double[m_dim];
      m_no_grid_points = new abs_type[m_dim];
      m_NN = new abs_type[m_dim];

      /* determine number grid points in each dimension */
      std::size_t no_l, no_u;
      int sign_l, sign_u;
      for(int index=0; index<m_dim; index++) {
        m_eta[index] = eta[index];
        /* ceil */
        try {
          /* get sign */
          sign_l = (lb[index] > 0) ? 1 : ((lb[index] < 0) ? -1 : 0);
          /* compute number of grid points from zero to lower bound */
          no_l=std::llround(std::abs(lb[index])/eta[index]+sign_l*0.5);
        } catch (...) {
          reset();
          std::ostringstream os;
          os << "\nscots::UniformGrid: something wrong in the division of " << lb[index] << " by " << eta[index] ;
          throw std::runtime_error(os.str().c_str());
        }
        /* floor */
        try {
          /* get sign */
          sign_u = (ub[index] > 0) ? 1 : ((ub[index] < 0) ? -1 : 0);
          /* compute number of grid points from zero to upper bound */
          no_u=std::llround(std::abs(ub[index])/eta[index]-sign_u*0.5);
        } catch (...) {
          reset();
          std::ostringstream os;
          os << "\nscots::UniformGrid: something wrong in the division of " << ub[index] << " by " << eta[index] ;
          throw std::runtime_error(os.str().c_str());
        }
        /* check if number of grid points in dimension index does not exceed max representable by abs_type  */
        if((sign_u*no_u-sign_l*no_l+1) > std::numeric_limits<abs_type>::max()) {
          reset();
          std::ostringstream os;
          throw std::runtime_error("\nscots::UniformGrid: number of grid points exceeds maximum value of abs_type (defined in TransitionSystem.hh).");
        }
        m_no_grid_points[index] = sign_u*no_u-sign_l*no_l+1;
        /* first grid point coordinates */
        m_first[index]= (double)sign_l*(double)no_l*eta[index];
      }
      /* compute m_NN */
      total_no_grid_points();
    } else {
      reset();
      throw std::runtime_error("\nscots::UniformGrid: grid dimension has to be greater than zero (using non-default constructor)");
    }
  }


  /** @brief compute the index associated with a grid point **/
  template<class grid_point_t>
  abs_type xtoi(const grid_point_t& x) const {
    abs_type idx = 0;
    double d_idx;
    double eta_h;

    for(int k=0; k<m_dim; k++) {
      d_idx = x[k]-m_first[k];
      eta_h = m_eta[k]/2.0;

      if ( d_idx <= -eta_h || d_idx >= m_no_grid_points[k]*m_eta[k]+eta_h ) {
        std::ostringstream os;
        os << "\nscots::UniformGrid: state ";
        for(int i=0; i<m_dim; i++) {
          os << x[i] << " ";
        }
        os << " is outside uniform grid.";
        throw std::runtime_error(os.str().c_str());
      }
      idx += static_cast<abs_type>((d_idx+eta_h )/m_eta[k])*m_NN[k];
    }
    return idx;
  }

  /** @brief compute the grid point associated with a index **/
  template<class grid_point_t>
  grid_point_t itox(abs_type idx) const {
    /* map index idx to grid point */
    grid_point_t x;
    abs_type num;
    for(abs_type k = m_dim-1; k > 0; k--) {
      num=idx/m_NN[k];
      idx=idx%m_NN[k];
      x[k]=m_first[k]+num*m_eta[k];
    }
    num=idx;
    x[0]=m_first[0]+num*m_eta[0];
    return x;
  }

  /** @brief creates console output with grid information **/
  void print_info(void) const {
    std::cout << "Distance of grid points (eta): ";
    for(int i=0; i<m_dim; i++) {
      std::cout << m_eta[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "Lower-left grid point: ";
    for(int i=0; i<m_dim; i++) {
      std::cout << m_first[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "Upper-right grid point: ";
    for(int i=0; i<m_dim; i++) {
      std::cout << m_first[i]+m_eta[i]*(m_no_grid_points[i]-1) << " ";
    }
    std::cout << std::endl;

    std::cout << "Number of grid points in each dimension: ";
    for(int i=0; i<m_dim; i++) {
        std::cout << m_no_grid_points[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "Number of grid points: "<< total_no_grid_points() << std::endl;
    std::cout << std::endl;
  }


  /* @brief function to write the grid into a file via a scots::FileWriter **/
//  bool addGridToFile(FileWriter&);            
  /* @brief function to read and reset(!) the grid from a file via a scots::FileReader **/
//  bool readFromGridFile(FileReader&, std::size_t); 

  /** @name get functions **/
  //@{
  int get_dim(void) const {
    return m_dim;
  }
  /* total number of grid points */
  abs_type size() const {
    return total_no_grid_points();
  }
  std::vector<double> get_eta(void) const {
    std::vector<double> eta;
    for(int i=0; i<m_dim; i++) {
      eta.push_back(m_eta[i]);
    }
    return eta;
  }
  std::vector<double> get_lower_left(void) const {
    std::vector<double> lower_left;
    for(int i=0; i<m_dim; i++) {
      lower_left.push_back(m_first[i]);
    }
    return lower_left;
  }
  std::vector<double> get_upper_right(void) const {
    std::vector<double> upper_right;
    for(int i=0; i<m_dim; i++) {
      upper_right.push_back(m_first[i]+m_eta[i]*(m_no_grid_points[i]-1));
    }
    return upper_right;
  }
  std::vector<abs_type> get_no_gp_per_dim() const {
    std::vector<abs_type> no_grid_points;
    for(int i=0; i<m_dim; i++) {
      no_grid_points.push_back(m_no_grid_points[i]);
    }
    return no_grid_points;
  }
  std::vector<abs_type> get_nn() const {
    std::vector<abs_type> NN;
    for(int i=0; i<m_dim; i++) {
      NN.push_back(m_NN[i]);
    }
    return NN;
  }
  //@}

private:
	/** @brief helper function to calculate the overall number of grid points and NN **/
  abs_type total_no_grid_points() {
    abs_type total=1;
    for(int i=0; i<m_dim; i++) {
      m_NN[i] = total;
      /* check overflow */
      if((total > std::numeric_limits<abs_type>::max()) || (total==std::numeric_limits<abs_type>::max() && m_no_grid_points[i]>1)) {
        reset();
        throw std::runtime_error("\nscots::UniformGrid: number of grid points exceeds maximum value of abs_type (defined in TransitionSystem.hh).");
      }
      total *= m_no_grid_points[i];
    }
    return total;
  }
	/** @brief helper for reseting the whole grid **/
  void reset();                            
	/** @brief setting everything to zero **/
  void null();                             
};


//bool UniformGrid::addGridToFile(FileWriter& writer) {
//  if(writer.open()) {
//    writer.add_TYPE(SCOTS_UG_TYPE);
//    writer.add_MEMBER(SCOTS_UG_DIM,m_dim);
//    writer.add_ARRAY(SCOTS_UG_GPPD,m_no_grid_points);
//    writer.add_ARRAY(SCOTS_UG_ETA,m_eta,m_dim);
//    writer.add_ARRAY(SCOTS_UG_FIRST,m_first,m_dim);
//
//    writer.close();
//    return true;
//  }
//  return false;
//}

//bool UniformGrid::readFromGridFile(FileReader& reader, std::size_t offset) {
//  reset();
//  if(!reader.open()) {
//    return false;
//  }
//  if(!reader.get_MEMBER(SCOTS_UG_DIM,m_dim,offset)) {
//    reset();
//    return false;
//  } else {
//   m_eta = new double[m_dim];
//   m_first = new double[m_dim];
//   m_no_grid_points = new abs_type[m_dim];
//   m_NN = new abs_type[m_dim];
//  }
//  if(!reader.get_ARRAY(SCOTS_UG_GPPD,m_no_grid_points,offset)) {
//    reset();
//    return false;
//  }
//  if(!reader.get_ARRAY(SCOTS_UG_ETA,m_eta,m_dim,offset)) {
//    reset();
//    return false;
//  }
//  if(!reader.get_ARRAY(SCOTS_UG_FIRST,m_first,m_dim,offset)) {
//    reset();
//    return false;
//  }
//
//  reader.close();
//  calculateTotalNumberOfGridPoints();
//  return true;
//}



void UniformGrid::reset() {
  m_dim = 0;

  delete[] m_eta;
  delete[] m_first;
  delete[] m_no_grid_points;
  delete[] m_NN;

  m_eta = nullptr;
  m_first = nullptr;
  m_no_grid_points = nullptr;
  m_NN = nullptr;
}

void UniformGrid::null() {
  m_dim = 0;
  m_eta = nullptr;
  m_first = nullptr;
  m_no_grid_points = nullptr;
  m_NN = nullptr;
}



///** 
// * @brief write atomic propositions to file
// */
//template<class F>
//bool writeAtomicPropositionsToFile(F& AP,abs_type numberOfCells,FileWriter& file,std::string name)
//{
//    if(!file.open())
//    {
//        return false;
//    }
//    scots::array<bool> table(numberOfCells);
//
//    for(size_t index = 0; index < numberOfCells; index++)
//    {
//        table[index] = AP(index);
//    }
//
//    return file.add_ARRAY<scots::array<bool> >(name,table,table.size());
//}

} /* close namespace */

#endif /* UNIFORMGRID_HH_ */
