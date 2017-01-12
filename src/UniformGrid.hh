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
#include <sstream>
#include <string>
#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <climits>
#include <set>

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
public:
  /* @cond  EXCLUDE from doxygen */
  /* default constructor */
  UniformGrid();                      
  /* destructor */
  ~UniformGrid();
  /* copy constructor */
  UniformGrid(const UniformGrid&);   
  /* move constructor */
  UniformGrid(UniformGrid&&);        
  /* copy assignment operator */
  UniformGrid& operator=(const UniformGrid&);
  /* move assignment operator */
  UniformGrid& operator=(UniformGrid&&);         
  /* @endcond */

  /* constructor */
  template<class grid_point_t>
  UniformGrid(const int,
              const grid_point_t&,
              const grid_point_t&,
              const grid_point_t&); 

  /** @brief compute the index associated with a grid point **/
  template<class grid_point_t>
  inline void xtoi(abs_type &, const grid_point_t&) const; 
  /** @brief compute the grid point associated with a index **/
  template<class grid_point_t>
  inline void itox(abs_type, grid_point_t &) const;  
  /** @brief creates console output with grid information **/
  void printInfo() const;              
  /* @brief function to write the grid into a file via a scots::FileWriter **/
//  bool addGridToFile(FileWriter&);            
  /* @brief function to read and reset(!) the grid from a file via a scots::FileReader **/
//  bool readFromGridFile(FileReader&, std::size_t); 

  /** @name get functions **/
  //@{
  inline int getDimension() const;
  inline abs_type getTotalNoGridPoints() const;
  inline std::vector<double> getEta() const;
  inline std::vector<double> getLowerLeftGridPoint() const;
  inline std::vector<double> getUpperRightGridPoint() const;
  inline std::vector<abs_type> getNoGridPointsPerDimension() const;
  inline std::vector<abs_type> getNN() const;
  //@}

private:
  /** @brief dimension of the Eucleadian space **/
  int m_dimension;                
  /** @brief m_dimension-dimensional vector containing the grid node distances **/
  double*      m_eta;                      
	/** @brief m_dimension-dimensional vector containing the real values of the first grid point **/
  double*      m_first_grid_point;        
	/** @brief scots::abs_type array[m_dimension] containing the number of grid points in each dimension **/
  abs_type*    m_no_grid_points;        
	/** @brief total number of grid points  **/
  abs_type     m_total_no_grid_points;  
	/** @brief array recursively defined by: m_NN[0]=1; m_NN[i]=m_NN[i-1}*no_grid_points[i-1]; **/
  abs_type*    m_NN;                       

	/** @brief helper for reseting the whole grid **/
  void reset();                            
	/** @brief setting everything to zero **/
  void null();                             
	/** @brief helper function to calculate the overall number of grid points and NN **/
  void calculateTotalNumberOfGridPoints(); 
};

UniformGrid::UniformGrid() {
  m_dimension = 0;
  m_total_no_grid_points = 0;
  m_eta = nullptr;
  m_first_grid_point = nullptr;
  m_no_grid_points = nullptr; 
  m_NN = nullptr;
}

UniformGrid::UniformGrid(const UniformGrid& other) : UniformGrid() {
  *this=other;
}

UniformGrid::UniformGrid(UniformGrid&& other) : UniformGrid() {
  *this=std::move(other);
}

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
UniformGrid::UniformGrid(const int dim, const grid_point_t& lb, const grid_point_t& ub, const grid_point_t& eta) : UniformGrid() {
  m_dimension = dim;
  if(m_dimension != 0) {
    /* check inut arguments */
    for(int index=0; index<dim; index++) {
      if(eta[index] <= 0) 
        throw std::runtime_error("\nscots::UniformGrid: eta must have positive entries.");
      if(lb[index] > ub[index]) 
        throw std::runtime_error("\nscots::UniformGrid: lower-left bound must be less than or equal to upper-right bound.");
    }
    m_eta = new double[m_dimension];
    m_first_grid_point = new double[m_dimension];
    m_no_grid_points = new abs_type[m_dimension];
    m_NN = new abs_type[m_dimension];

    /* determine number grid points in each dimension */
    std::size_t no_l, no_u;
    int sign_l, sign_u;
    for(int index=0; index<m_dimension; index++) {
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
      m_first_grid_point[index]= (double)sign_l*(double)no_l*eta[index];
    }
    /* total number of total grid points */
    calculateTotalNumberOfGridPoints();
  } else {
    reset();
    throw std::runtime_error("\nscots::UniformGrid: grid dimension has to be greater than zero (using non-default constructor)");
  }
}

UniformGrid::~UniformGrid() {
  delete[] m_eta;
  delete[] m_first_grid_point;
  delete[] m_no_grid_points;
  delete[] m_NN;
}

UniformGrid &UniformGrid::operator=(const UniformGrid &other) {
	if(this==&other)
		return *this;
	reset();
  m_dimension=other.m_dimension;
  m_total_no_grid_points=other.m_total_no_grid_points;
  if(m_dimension != 0) {
    m_eta = new double[m_dimension];
    m_first_grid_point = new double[m_dimension];
    m_no_grid_points = new abs_type[m_dimension];
    m_NN = new abs_type[m_dimension];
    for(int index=0; index<m_dimension; index++) {
      m_eta[index] = other.m_eta[index];
      m_first_grid_point[index]  = other.m_first_grid_point[index];
      m_no_grid_points[index]  = other.m_no_grid_points[index];
      m_NN[index]  = other.m_NN[index];
    }
  } 
  return *this;
}

UniformGrid& UniformGrid::operator=(UniformGrid&& other) {
	reset();

  m_dimension=other.m_dimension;
	m_eta=other.m_eta;
	m_first_grid_point=other.m_first_grid_point;
	m_no_grid_points=other.m_no_grid_points;
	m_NN=other.m_NN;
	m_total_no_grid_points=other.m_total_no_grid_points;

  other.null();
	return *this;
}

template<class grid_point_t>
inline void UniformGrid::xtoi(abs_type &idx, const grid_point_t& x) const {
  idx = 0;
  double d_idx;
  double eta_h;

  for(int k=0; k<m_dimension; k++) {
    d_idx = x[k]-m_first_grid_point[k];
    eta_h = m_eta[k]/2.0;

    if ( d_idx <= -eta_h || d_idx >= m_no_grid_points[k]*m_eta[k]+eta_h ) {
      std::ostringstream os;
      os << "\nscots::UniformGrid: x is outside uniform grid." << x[k] ;
      throw std::runtime_error(os.str().c_str());
    }
    idx += static_cast<abs_type>((d_idx+eta_h )/m_eta[k])*m_NN[k];
  }
}

template<class grid_point_t>
inline void UniformGrid::itox(abs_type idx, grid_point_t &x) const {
  if(idx >= m_total_no_grid_points) {
    std::ostringstream os;
    os << "\nscots::UniformGrid: idx larger than number of grid points.";
    throw std::runtime_error(os.str().c_str());
  }
  /* map index idx to grid point */
  abs_type num;
  for(abs_type k = m_dimension-1; k > 0; k--) {
    num=idx/m_NN[k];
    idx=idx%m_NN[k];
    x[k]=m_first_grid_point[k]+num*m_eta[k];
  }
  num=idx;
  x[0]=m_first_grid_point[0]+num*m_eta[0];
}

void UniformGrid::printInfo(void) const {
  std::cout << "Distance of grid points (eta): ";
  for(int i=0; i<m_dimension; i++) {
    std::cout << m_eta[i] << " ";
  }
  std::cout << std::endl;

  std::cout << "Lower-left grid point: ";
  for(int i=0; i<m_dimension; i++) {
    std::cout << m_first_grid_point[i] << " ";
  }
  std::cout << std::endl;

  std::cout << "Upper-right grid point: ";
  for(int i=0; i<m_dimension; i++) {
    std::cout << m_first_grid_point[i]+m_eta[i]*(m_no_grid_points[i]-1) << " ";
  }
  std::cout << std::endl;

  std::cout << "Number of grid points in each dimension: ";
  for(int i=0; i<m_dimension; i++) {
      std::cout << m_no_grid_points[i] << " ";
  }
  std::cout << std::endl;

  std::cout << "Number of grid points: "<< m_total_no_grid_points << std::endl;
  std::cout << std::endl;
}

//bool UniformGrid::addGridToFile(FileWriter& writer) {
//  if(writer.open()) {
//    writer.add_TYPE(SCOTS_UG_TYPE);
//    writer.add_MEMBER(SCOTS_UG_DIM,m_dimension);
//    writer.add_ARRAY(SCOTS_UG_GPPD,m_no_grid_points);
//    writer.add_ARRAY(SCOTS_UG_ETA,m_eta,m_dimension);
//    writer.add_ARRAY(SCOTS_UG_FIRST,m_first_grid_point,m_dimension);
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
//  if(!reader.get_MEMBER(SCOTS_UG_DIM,m_dimension,offset)) {
//    reset();
//    return false;
//  } else {
//   m_eta = new double[m_dimension];
//   m_first_grid_point = new double[m_dimension];
//   m_no_grid_points = new abs_type[m_dimension];
//   m_NN = new abs_type[m_dimension];
//  }
//  if(!reader.get_ARRAY(SCOTS_UG_GPPD,m_no_grid_points,offset)) {
//    reset();
//    return false;
//  }
//  if(!reader.get_ARRAY(SCOTS_UG_ETA,m_eta,m_dimension,offset)) {
//    reset();
//    return false;
//  }
//  if(!reader.get_ARRAY(SCOTS_UG_FIRST,m_first_grid_point,m_dimension,offset)) {
//    reset();
//    return false;
//  }
//
//  reader.close();
//  calculateTotalNumberOfGridPoints();
//  return true;
//}


inline int UniformGrid::getDimension(void) const {
  return m_dimension;
}

inline abs_type UniformGrid::getTotalNoGridPoints() const {
  return m_total_no_grid_points;
}

inline std::vector<double> UniformGrid::getEta(void) const {
  std::vector<double> eta;
  for(int i=0; i<m_dimension; i++) {
    eta.push_back(m_eta[i]);
  }
  return eta;
}

inline std::vector<double> UniformGrid::getLowerLeftGridPoint(void) const {
  std::vector<double> lower_left;
  for(int i=0; i<m_dimension; i++) {
    lower_left.push_back(m_first_grid_point[i]);
  }
  return lower_left;
}

inline std::vector<double> UniformGrid::getUpperRightGridPoint(void) const {
  std::vector<double> upper_right;
  for(int i=0; i<m_dimension; i++) {
    upper_right.push_back(m_first_grid_point[i]+m_eta[i]*(m_no_grid_points[i]-1));
  }
  return upper_right;
}
inline std::vector<abs_type> UniformGrid::getNoGridPointsPerDimension() const {
  std::vector<abs_type> no_grid_points;
  for(int i=0; i<m_dimension; i++) {
    no_grid_points.push_back(m_no_grid_points[i]);
  }
  return no_grid_points;
}

inline std::vector<abs_type> UniformGrid::getNN() const {
  std::vector<abs_type> NN;
  for(int i=0; i<m_dimension; i++) {
    NN.push_back(m_NN[i]);
  }
  return NN;
}

void UniformGrid::reset() {
  m_dimension = 0;
  m_total_no_grid_points = 0;

  delete[] m_eta;
  delete[] m_first_grid_point;
  delete[] m_no_grid_points;
  delete[] m_NN;

  m_eta = nullptr;
  m_first_grid_point = nullptr;
  m_no_grid_points = nullptr;
  m_NN = nullptr;
}

void UniformGrid::null() {
  m_dimension = 0;
  m_total_no_grid_points = 0;
  m_eta = nullptr;
  m_first_grid_point = nullptr;
  m_no_grid_points = nullptr;
  m_NN = nullptr;
}

void UniformGrid::calculateTotalNumberOfGridPoints() {
  std::size_t total=1;
  for(int i=0; i<m_dimension; i++) {
    m_NN[i] = total;
    /* check overflow */
    if((total > std::numeric_limits<abs_type>::max()) || (total==std::numeric_limits<abs_type>::max() && m_no_grid_points[i]>1)) {
      reset();
      throw std::runtime_error("\nscots::UniformGrid: number of grid points exceeds maximum value of abs_type (defined in TransitionSystem.hh).");
    }
    total *= m_no_grid_points[i];
  }
  m_total_no_grid_points = total;
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
