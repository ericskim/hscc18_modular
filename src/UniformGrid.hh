/*
 * UniformGrid.hh
 *
 *  created on: 26.10.2015
 *      author: Matthias Rungger
 *              Frederik Kunik
 */

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

//#include "TransitionSystem.hh"
//#include "FileHandler.hh"

using abs_type=uint32_t;

#define SCOTS_UG_TYPE   "UNIFORMGRID"
#define SCOTS_UG_DIM    "DIM"
#define SCOTS_UG_FIRST  "FIRST"
#define SCOTS_UG_ETA    "ETA"
#define SCOTS_UG_GPPD   "GPPD"


namespace scots {

/*!
 * \brief The UniformGrid class stores information of a uniform grid \n
 *
 * Properties: \n
 * - the grid points are distributied uniformly in each dimension \n
 * - the domain of the unfiorm grid is defined by a hyper interval \n
 * - each grid pont is associated with a cell, i.e. a hyper rectangle with \n
 *   radius eta/2 centered at the grid point \n
 *
 * Grid point alignment: \n
 * - the origin is a grid point (not necessarily contained in the set) \n
 * - the distance of the grid points in each dimension i is defined by eta[i] \n
 *
 * See \n
 * - the manual in 
 * - http://arxiv.org/abs/1503.03715 for theoretical background \n
 *
 *
 */
class UniformGrid {
public:
  UniformGrid();                      //!< default contructor
  UniformGrid(const UniformGrid&);    //!< copy contructor
  UniformGrid(UniformGrid&&);         //!< move contructor
  template<class grid_point_t>
  UniformGrid(const unsigned int,
              const grid_point_t&,
              const grid_point_t&,
              const grid_point_t&); //!< non-default constructor

  ~UniformGrid();

  UniformGrid& operator=(const UniformGrid&);    //!< copy asignment operator
  UniformGrid& operator=(UniformGrid&&);         //!< move asignment operator

  template<class grid_point_t>
  inline void xtoi(abs_type &, const grid_point_t&) const; //!< compute the index associated with a grid point
  template<class grid_point_t>
  inline void itox(abs_type, grid_point_t &) const;       //!< compute the grid point associated with an index

  void printInfo() const;              //!< creates console output with grid information
//  bool addGridToFile(FileWriter&);            //!< function to write the grid into a file via a scots::FileWriter
//  bool readFromGridFile(FileReader&, std::size_t); //!< function to read and reset(!) the grid from a file via a scots::FileReader

  /* get functions */
  inline unsigned int getDimension() const;
  inline abs_type getTotalNoOfGridPoints() const;
  inline std::vector<double> getEta() const;
  inline std::vector<double> getFirstGridPoint() const;
  inline std::vector<abs_type> getNoOfGridPointsPerDimension() const;
  inline std::vector<abs_type> getNN() const;

private:
  unsigned int m_dimension;                //!< dimension of the real space
  double*      m_eta;                      //!< m_dimension-dimensional vector containing the grid node distances
  double*      m_first_grid_point;        //!< m_dimension-dimensional vector containing the real values of the first grid point
  abs_type*    m_no_of_grid_points;        //!< scots::abs_type array[m_dimension] containing the number of grid points in each dimension */
  abs_type     m_total_no_of_grid_points;  //!< total number of grid points */
  abs_type*    m_NN;                       //!< ToDo

  void reset();                            //!< helper for reseting the whole grid
  void calculateTotalNumberOfGridPoints(); //!< helper function to calculate the overall number of grid points and NN
};

UniformGrid::UniformGrid() {
  m_dimension = 0;
  m_total_no_of_grid_points = 0;
  m_eta = NULL;
  m_first_grid_point = NULL;
  m_no_of_grid_points = NULL; 
  m_NN = NULL;
}

UniformGrid::UniformGrid(const UniformGrid& other) : UniformGrid() {
  *this=other;
}

UniformGrid::UniformGrid(UniformGrid&& other) : UniformGrid() {
  *this=std::move(other);
}
/*!
 * \brief UniformGrid::UniformGrid
 *  provide uniform grid parameters and domain defining hyper interval \n
 * \param dim   - dimension of the real space
 * \param lb    - lower-left corner of the hyper-interval confining the uniform grid
 * \param ub    - upper-right corner of the hyper-interval confining the uniform grid
 * \param eta   - grid point distances
 */
template<class grid_point_t>
UniformGrid::UniformGrid(const unsigned int dim, const grid_point_t& lb, const grid_point_t& ub, const grid_point_t& eta) : UniformGrid::UniformGrid() {
  m_dimension = dim;
  if(m_dimension != 0) {
    /* check inut arguments */
    for(std::size_t index=0; index<dim; index++) {
      if(eta[index] <= 0) 
        throw std::runtime_error("\nscots::UniformGrid: eta must have positive entries.");
      if(lb[index] > ub[index]) 
        throw std::runtime_error("\nscots::UniformGrid: lower-left bound must be less than or equal to upper-right bound.");
    }
    m_eta = new double[m_dimension];
    m_first_grid_point = new double[m_dimension];
    m_no_of_grid_points = new abs_type[m_dimension];
    m_NN = new abs_type[m_dimension];

    /* determine number grid points in each dimension */
    std::size_t no_l, no_u;
    int sign_l, sign_u;
    for(std::size_t index=0; index<m_dimension; index++) {
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
      m_no_of_grid_points[index] = sign_u*no_u-sign_l*no_l+1;
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
  delete[] m_no_of_grid_points;
  delete[] m_NN;
}

UniformGrid &UniformGrid::operator=(const UniformGrid &other) {
	if(this==&other)
		return *this;
	reset();
  m_dimension=other.m_dimension;
  if(m_dimension != 0) {
    m_eta = new double[m_dimension];
    m_first_grid_point = new double[m_dimension];
    m_no_of_grid_points = new abs_type[m_dimension];
    m_NN = new abs_type[m_dimension];
    for(std::size_t index=0; index<m_dimension; index++) {
      m_eta[index] = other.m_eta[index];
      m_first_grid_point[index]  = other.m_first_grid_point[index];
      m_no_of_grid_points[index]  = other.m_no_of_grid_points[index];
      m_NN[index]  = other.m_NN[index];
    }
  } 
  return *this;
}


UniformGrid& UniformGrid::operator=(UniformGrid&& other) {
	if(this==&other)
		return *this;
	reset();

  m_dimension=other.m_dimension;
	m_eta=other.m_eta;
	m_first_grid_point=other.m_first_grid_point;
	m_no_of_grid_points=other.m_no_of_grid_points;
	m_NN=other.m_NN;
	m_total_no_of_grid_points=other.m_total_no_of_grid_points;

	other.reset();

	return *this;
}

template<class grid_point_t>
inline void UniformGrid::xtoi(abs_type &idx, const grid_point_t& x) const {
  idx = 0;
  double d_idx;
  double eta_h;

  for(std::size_t k=0; k<m_dimension; k++) {
    d_idx = x[k]-m_first_grid_point[k];
    eta_h = m_eta[k]/2.0;

    if ( d_idx <= -eta_h || d_idx >= m_no_of_grid_points[k]*m_eta[k]+eta_h ) {
      std::ostringstream os;
      os << "\nscots::UniformGrid: x is outside uniform grid." << x[k] ;
      throw std::runtime_error(os.str().c_str());
    }
    idx += static_cast<abs_type>((d_idx+eta_h )/m_eta[k])*m_NN[k];
  }
}

template<class grid_point_t>
inline void UniformGrid::itox(abs_type idx, grid_point_t &x) const {
  if(idx >= m_total_no_of_grid_points) {
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
  std::cout << "Grid node distance (eta) in each dimension: ";
  for(std::size_t i=0; i<m_dimension; i++) {
    std::cout << m_eta[i] << " ";
  }
  std::cout << std::endl;

  std::cout << "First grid point: ";
  for(std::size_t i=0; i<m_dimension; i++) {
    std::cout << m_first_grid_point[i] << " ";
  }
  std::cout << std::endl;

  std::cout << "Number of grid points in each dimension: ";
  for(std::size_t i=0; i<m_dimension; i++) {
      std::cout << m_no_of_grid_points[i] << " ";
  }
  std::cout << std::endl;

  std::cout << "Number of grid points: "<< m_total_no_of_grid_points << std::endl;
}

//bool UniformGrid::addGridToFile(FileWriter& writer) {
//  if(writer.open()) {
//    writer.add_TYPE(SCOTS_UG_TYPE);
//    writer.add_MEMBER(SCOTS_UG_DIM,m_dimension);
//    writer.add_ARRAY(SCOTS_UG_GPPD,m_no_of_grid_points);
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
//   m_no_of_grid_points = new abs_type[m_dimension];
//   m_NN = new abs_type[m_dimension];
//  }
//  if(!reader.get_ARRAY(SCOTS_UG_GPPD,m_no_of_grid_points,offset)) {
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


inline abs_type UniformGrid::getDimension(void) const {
  return m_dimension;
}

inline abs_type UniformGrid::getTotalNoOfGridPoints() const {
  return m_total_no_of_grid_points;
}

inline std::vector<double> UniformGrid::getEta(void) const {
  std::vector<double> eta(m_dimension);
  for(std::size_t index=0; index<m_dimension; index++) {
    eta[index]=m_eta[index];
  }
  return std::move(eta);
}

inline std::vector<double> UniformGrid::getFirstGridPoint(void) const {
  std::vector<double> first(m_dimension);
  for(std::size_t index=0; index<m_dimension; index++) {
    first[index]=m_first_grid_point[index];
  }
  return std::move(first);
}

inline std::vector<abs_type> UniformGrid::getNoOfGridPointsPerDimension() const {
  std::vector<abs_type> no_of_grid_points(m_dimension);
  for(std::size_t index=0; index<m_dimension; index++) {
    no_of_grid_points.push_back(m_no_of_grid_points[index]);
  }
  return std::move(no_of_grid_points);
}

inline std::vector<abs_type> UniformGrid::getNN() const {
  std::vector<abs_type> NN(m_dimension);
  for(std::size_t index=0; index<m_dimension; index++) {
    NN.push_back(m_NN[index]);
  }
  return std::move(NN);
}

void UniformGrid::reset() {
  m_dimension = 0;
  m_total_no_of_grid_points = 0;

  delete[] m_eta;
  delete[] m_first_grid_point;
  delete[] m_no_of_grid_points;
  delete[] m_NN;

  m_eta = NULL;
  m_first_grid_point = NULL;
  m_no_of_grid_points = NULL;
  m_NN = NULL;
}

void UniformGrid::calculateTotalNumberOfGridPoints() {
  std::size_t total=1;
  for(std::size_t index=0; index<m_dimension; index++) {
    m_NN[index] = total;
    /* check overflow */
    if((total > std::numeric_limits<abs_type>::max()) || (total==std::numeric_limits<abs_type>::max() && m_no_of_grid_points[index]>1)) {
      reset();
      throw std::runtime_error("\nscots::UniformGrid: number of grid points exceeds maximum value of abs_type (defined in TransitionSystem.hh).");
    }
    total *= m_no_of_grid_points[index];
  }
  m_total_no_of_grid_points = total;
}

} /* close namespace */

#endif /* UNIFORMGRID_HH_ */
